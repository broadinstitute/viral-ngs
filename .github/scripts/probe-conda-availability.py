#!/usr/bin/env python3
"""Probe conda channels for the availability of CVE fix versions.

Written for the Claude triage / auto-fix pipeline in `.github/workflows/container-scan.yml`.
Those agents are asked whether a version-floor bump in `docker/requirements/*.txt` can
close a CVE, which requires knowing:

  (a) whether the fixed version exists on conda-forge or bioconda at all (conda-forge
      often lags PyPI, so a Trivy `FixedVersion` is not proof of availability), and
  (b) whether it is built for every architecture we publish -- `linux-64` AND
      `linux-aarch64`, or `noarch`.

Historically the agents had no tool that could answer either question, so they skipped
every floor bump as "unverifiable" -- see the SKIPPED.md decision logs from runs
29883237168 (Pillow) and 30898278357 (cryptography). This script answers both questions
deterministically, before the agent runs, and writes the result to a JSON file the agent
reads. Keeping the lookup here rather than leaving it to the agent's own web access means
the answer is identical run-to-run, costs the agent no turns, and is archived as a
workflow artifact for audit.

It also reports the reverse direction: which *other* packages pinned in
`docker/requirements/*.txt` carry a dependency constraint that would exclude the fix
version. That is the trap CVE-2026-69247 hit -- `pyopenssl` 26.3.0 pins
`cryptography >=49.0.0,<50`, so bumping `cryptography` alone requires `pyopenssl` 26.4.0
(which widened the cap to `<51`) to come along.

Usage:
    probe-conda-availability.py \
        --trivy-json trivy-results.json \
        --cve-ids "CVE-2026-69247 CVE-2025-47273" \
        --requirements-dir docker/requirements \
        --out .cve-fix-context/conda-availability.json

Exit status is 0 even when lookups fail: partial data plus a populated `probe_errors`
list is more useful to the agent than a failed workflow step. The agent is instructed to
treat a package with no usable data as unverified -- i.e. SKIP -- which is the safe
direction.
"""

import argparse
import concurrent.futures
import glob
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request

ANACONDA_API = "https://api.anaconda.org/package/{channel}/{name}"
CHANNELS = ("conda-forge", "bioconda")

# Architectures the published images are built for; see the build matrix in docker.yml.
REQUIRED_SUBDIRS = ("linux-64", "linux-aarch64")
# The conda env in the images is Python 3.12 (see docker/Dockerfile.baseimage).
PYTHON_TAG = "py312"

# How many versions of a constraining package to walk when looking for the oldest one
# that admits the fix version. Bounded so a package with a decade of releases doesn't
# turn into a pathological loop.
MAX_VERSIONS_SCANNED = 60

# The reverse scan is a bonus, not the critical path, and some payloads are large
# (google-cloud-sdk is ~7 MB). Fetch candidates concurrently and give the whole scan a
# wall-clock budget; if it runs out we report what we have plus a probe_error saying the
# scan was truncated, rather than silently returning partial data.
REVERSE_SCAN_WORKERS = 8
REVERSE_SCAN_BUDGET_SECONDS = 90
# Reverse-scan fetches get one shot with a short timeout: a package whose payload is too
# big to read in time (awscli lists thousands of releases) will time out again on retry,
# and burning the budget on it starves the rest of the scan. Unread packages are reported
# in `probe_errors` rather than silently dropped.
REVERSE_SCAN_TIMEOUT = 12

USER_AGENT = "viral-ngs-container-scan/1.0 (+https://github.com/broadinstitute/viral-ngs)"


# --------------------------------------------------------------------------------------
# Version handling
#
# Prefer `packaging` when the runner has it; fall back to a local comparator otherwise so
# this script has no dependencies outside the stdlib. Conda-style versions such as
# `2026.6.3` (rpds-py) must parse under either path.
# --------------------------------------------------------------------------------------

try:
    from packaging.version import InvalidVersion, Version

    def version_key(raw):
        try:
            return (0, Version(raw))
        except InvalidVersion:
            return (1, _fallback_key(raw))

except ImportError:  # pragma: no cover - depends on runner image

    def version_key(raw):
        return (1, _fallback_key(raw))


def _fallback_key(raw):
    """Split a version into comparable (kind, value) chunks.

    Numeric chunks sort as integers so 10 > 9; alphabetic chunks sort after numeric ones
    at the same position. Good enough for "is this release >= that release", which is all
    we need -- we are not trying to model conda's full pre-release ordering.
    """
    key = []
    for chunk in re.split(r"[._\-+]", str(raw)):
        for match in re.finditer(r"\d+|[A-Za-z]+", chunk):
            token = match.group()
            key.append((0, int(token), "") if token.isdigit() else (1, 0, token.lower()))
    return tuple(key)


def version_ge(a, b):
    return version_key(a) >= version_key(b)


def version_gt(a, b):
    return version_key(a) > version_key(b)


def sorted_versions(versions):
    return sorted(versions, key=version_key)


# --------------------------------------------------------------------------------------
# Conda constraint evaluation
# --------------------------------------------------------------------------------------

_OP_RE = re.compile(r"^(>=|<=|==|!=|>|<|=)?\s*(.+)$")


def constraint_allows(constraint, target):
    """Does a conda version constraint admit `target`?

    Handles the comma-separated conjunctions and `|` disjunctions that appear in conda
    `depends` strings (e.g. `>=49.0.0,<51`). Returns None when the constraint cannot be
    parsed, so callers can report uncertainty instead of guessing.
    """
    if not constraint or constraint.strip() in ("*", ""):
        return True

    # `|` is a disjunction in conda match specs; any satisfied branch admits the target.
    if "|" in constraint:
        results = [constraint_allows(part, target) for part in constraint.split("|")]
        if any(r is True for r in results):
            return True
        return None if any(r is None for r in results) else False

    for clause in constraint.split(","):
        clause = clause.strip()
        if not clause:
            continue
        match = _OP_RE.match(clause)
        if not match:
            return None
        op, bound = match.group(1) or "==", match.group(2).strip()

        # Prefix wildcards: `1.2.*` means "starts with 1.2.".
        if bound.endswith(".*") or bound.endswith("*"):
            prefix = bound.rstrip("*").rstrip(".")
            hit = str(target) == prefix or str(target).startswith(prefix + ".")
            if op in ("==", "="):
                if not hit:
                    return False
                continue
            if op == "!=":
                if hit:
                    return False
                continue
            return None

        try:
            if op == ">=":
                ok = version_ge(target, bound)
            elif op == ">":
                ok = version_gt(target, bound)
            elif op == "<=":
                ok = version_ge(bound, target)
            elif op == "<":
                ok = version_gt(bound, target)
            elif op in ("==", "="):
                ok = version_key(target) == version_key(bound)
            elif op == "!=":
                ok = version_key(target) != version_key(bound)
            else:
                return None
        except Exception:
            return None
        if not ok:
            return False
    return True


def constraint_floor(constraint):
    """Lowest version a constraint admits, as far as we can tell from `>=`/`>`/`==`."""
    if not constraint:
        return None
    best = None
    for clause in constraint.split(","):
        match = _OP_RE.match(clause.strip())
        if not match:
            continue
        op, bound = match.group(1) or "==", match.group(2).strip()
        if op in (">=", ">", "==", "=") and not bound.endswith("*"):
            if best is None or version_gt(bound, best):
                best = bound
    return best


# --------------------------------------------------------------------------------------
# Requirements files
# --------------------------------------------------------------------------------------

_REQ_LINE_RE = re.compile(r"^([A-Za-z0-9_.\-]+)\s*(.*)$")


def read_requirements(requirements_dir):
    """Map conda package name -> {file, line, spec, constraint} from docker/requirements."""
    pins = {}
    for path in sorted(glob.glob(os.path.join(requirements_dir, "*.txt"))):
        try:
            with open(path, encoding="utf-8") as handle:
                lines = handle.readlines()
        except OSError:
            continue
        for lineno, raw in enumerate(lines, start=1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            match = _REQ_LINE_RE.match(line)
            if not match:
                continue
            name, constraint = match.group(1), match.group(2).strip()
            pins.setdefault(
                name.lower(),
                {
                    "package": name,
                    "file": path,
                    "line": lineno,
                    "spec": line,
                    "constraint": constraint or None,
                },
            )
    return pins


# --------------------------------------------------------------------------------------
# anaconda.org lookups
# --------------------------------------------------------------------------------------


class ChannelClient:
    """Tiny cached client for the anaconda.org package API."""

    def __init__(self, errors, timeout=30, retries=3):
        self._cache = {}
        self._errors = errors
        self._timeout = timeout
        self._retries = retries
        self._lock = threading.Lock()

    def cached(self, channel, name):
        with self._lock:
            return self._cache.get((channel, name.lower()), False)

    def fetch(self, channel, name, timeout=None, retries=None):
        key = (channel, name.lower())
        with self._lock:
            if key in self._cache:
                return self._cache[key]

        timeout = self._timeout if timeout is None else timeout
        retries = self._retries if retries is None else retries
        url = ANACONDA_API.format(channel=channel, name=name)
        result = None
        for attempt in range(1, retries + 1):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    result = json.load(response)
                break
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    break  # Package simply isn't in this channel; not an error.
                if exc.code in (429, 500, 502, 503, 504) and attempt < retries:
                    time.sleep(2 * attempt)
                    continue
                self._errors.append(
                    "HTTP %s fetching %s/%s" % (exc.code, channel, name)
                )
                break
            except Exception as exc:  # network error, malformed JSON, ...
                if attempt < retries:
                    time.sleep(2 * attempt)
                    continue
                self._errors.append(
                    "%s fetching %s/%s: %s -- this package was NOT checked; treat any "
                    "claim about it as unverified"
                    % (type(exc).__name__, channel, name, exc)
                )
                break

        with self._lock:
            self._cache[key] = result
        return result


def files_by_version(payload):
    """Group a package payload's files by version -> list of {subdir, build, depends}."""
    grouped = {}
    for entry in payload.get("files") or []:
        attrs = entry.get("attrs") or {}
        version = entry.get("version")
        if not version:
            continue
        grouped.setdefault(version, []).append(
            {
                "subdir": attrs.get("subdir") or entry.get("subdir"),
                "build": attrs.get("build") or "",
                "depends": attrs.get("depends") or [],
                "upload_time": entry.get("upload_time"),
            }
        )
    return grouped


def describe_version(version, builds):
    subdirs = sorted({b["subdir"] for b in builds if b["subdir"]})
    noarch = "noarch" in subdirs
    py_tagged = {
        b["subdir"]
        for b in builds
        if b["subdir"] and b["build"].startswith(PYTHON_TAG)
    }
    uploads = sorted(b["upload_time"] for b in builds if b["upload_time"])
    return {
        "version": version,
        "subdirs": subdirs,
        "noarch": noarch,
        "linux_64": noarch or "linux-64" in subdirs,
        "linux_aarch64": noarch or "linux-aarch64" in subdirs,
        "py312_linux_64": noarch or "linux-64" in py_tagged,
        "py312_linux_aarch64": noarch or "linux-aarch64" in py_tagged,
        "arch_ok": noarch or all(s in subdirs for s in REQUIRED_SUBDIRS),
        "first_upload": uploads[0] if uploads else None,
    }


def probe_channel(client, channel, pkg_name, threshold):
    """Availability of `pkg_name` at versions >= `threshold` in one channel."""
    payload = client.fetch(channel, pkg_name)
    name_lookup = "exact"
    if payload is None and pkg_name != pkg_name.lower():
        payload = client.fetch(channel, pkg_name.lower())
        name_lookup = "lowercased"
    if payload is None:
        return {"name_lookup": "not_found"}

    grouped = files_by_version(payload)
    candidates = []
    for version in sorted_versions(grouped):
        if threshold and not version_ge(version, threshold):
            continue
        candidates.append(describe_version(version, grouped[version]))

    usable = [c for c in candidates if c["arch_ok"]]
    return {
        "name_lookup": name_lookup,
        "resolved_name": payload.get("name", pkg_name),
        "latest_version": payload.get("latest_version"),
        "candidates": candidates,
        "recommended_floor": usable[0]["version"] if usable else None,
        "arch_ok": bool(usable),
    }


def reverse_scan_candidates(pins, target_name):
    """Which pinned packages are worth checking for a constraint on `target_name`.

    Scanning all ~88 pinned packages costs a lot of large payloads for little return, so
    prefer the requirements file that pins the target -- a cap on a Python security lib
    lives alongside it in baseimage.txt, which is exactly where the `pyopenssl` cap on
    `cryptography` lives. Fall back to everything when the target isn't pinned at all.
    """
    target_lower = target_name.lower()
    pin = pins.get(target_lower)
    if pin:
        same_file = {
            key: value
            for key, value in pins.items()
            if value["file"] == pin["file"] and key != target_lower
        }
        if same_file:
            return same_file, pin["file"]
    return {k: v for k, v in pins.items() if k != target_lower}, None


def find_reverse_constraints(client, candidates, target_name, target_version, errors):
    """Packages pinned in docker/requirements that constrain `target_name`.

    This is the check that would have caught `pyopenssl <50` blocking `cryptography`
    50.0.0 without a human noticing it.
    """
    findings = []
    target_lower = target_name.lower()
    deadline = time.monotonic() + REVERSE_SCAN_BUDGET_SECONDS

    # Warm the cache concurrently: conda-forge first, then bioconda only for the misses.
    for channel in CHANNELS:
        pending = [
            pin["package"]
            for key, pin in sorted(candidates.items())
            if client.cached(channel, pin["package"]) is False
            and not any(
                client.cached(prior, pin["package"]) for prior in CHANNELS[: CHANNELS.index(channel)]
            )
        ]
        if not pending or time.monotonic() > deadline:
            continue
        with concurrent.futures.ThreadPoolExecutor(REVERSE_SCAN_WORKERS) as pool:
            futures = {
                pool.submit(
                    client.fetch,
                    channel,
                    name,
                    timeout=REVERSE_SCAN_TIMEOUT,
                    retries=1,
                ): name
                for name in pending
            }
            for future in concurrent.futures.as_completed(
                futures, timeout=max(1.0, deadline - time.monotonic())
            ):
                try:
                    future.result()
                except concurrent.futures.TimeoutError:
                    break
                except Exception:
                    pass  # ChannelClient already recorded the error.

    unscanned = []
    for key, pin in sorted(candidates.items()):
        for channel in CHANNELS:
            payload = client.cached(channel, pin["package"])
            if payload is False:
                # Never fetched (budget ran out before this one).
                if pin["package"] not in unscanned:
                    unscanned.append(pin["package"])
                continue
            if payload is None:
                continue

            grouped = files_by_version(payload)
            latest = payload.get("latest_version")

            def constraint_for(version):
                for build in grouped.get(version, []):
                    for dep in build["depends"]:
                        parts = dep.split(None, 1)
                        if parts and parts[0].lower() == target_lower:
                            return dep, (parts[1].strip() if len(parts) > 1 else "")
                return None, None

            latest_dep, latest_constraint = constraint_for(latest)
            if latest_dep is None:
                break  # This package doesn't depend on the target at all.

            latest_allows = constraint_allows(latest_constraint, target_version)

            # Walk newest-first to find the oldest version that still admits the target,
            # so we can say "raise the floor to exactly this".
            ordered = list(reversed(sorted_versions(grouped)))[:MAX_VERSIONS_SCANNED]
            min_allowing = None
            for version in ordered:
                _, constraint = constraint_for(version)
                if constraint is None:
                    continue
                if constraint_allows(constraint, target_version) is True:
                    min_allowing = version
                else:
                    break

            repo_floor = constraint_floor(pin["constraint"])
            needs_raise = bool(
                min_allowing and repo_floor and version_gt(min_allowing, repo_floor)
            ) or bool(min_allowing and not repo_floor)

            if latest_allows is None:
                errors.append(
                    "could not parse %s constraint %r on %s"
                    % (pin["package"], latest_constraint, target_name)
                )

            findings.append(
                {
                    "package": pin["package"],
                    "channel": channel,
                    "latest_version": latest,
                    "constraint_at_latest": latest_dep,
                    "latest_allows_target": latest_allows,
                    "min_version_allowing_target": min_allowing,
                    "repo_file": pin["file"],
                    "repo_line": pin["line"],
                    "repo_spec": pin["spec"],
                    "repo_floor": repo_floor,
                    "floor_needs_raise": needs_raise,
                }
            )
            break  # Found it in this channel; don't double-report from the other.

    if unscanned:
        errors.append(
            "reverse dependency scan for %s was truncated after %ss; %d package(s) not "
            "checked for a constraint on %s: %s"
            % (
                target_name,
                REVERSE_SCAN_BUDGET_SECONDS,
                len(unscanned),
                target_name,
                ", ".join(unscanned),
            )
        )

    return findings


# --------------------------------------------------------------------------------------
# Trivy input
# --------------------------------------------------------------------------------------


def classify_fix_shape(trivy_class, trivy_type):
    """What kind of fix a finding needs -- not every CVE is a conda floor bump.

    Without this, an Ubuntu `curl` finding whose fix version is `8.14.1-2ubuntu1.4` would
    be matched against conda-forge's `curl` 8.16.0 and reported as "floor bump verified",
    which is wrong: an apt package is not fixed by pinning a conda package of the same
    name. Probing conda for os-pkgs findings is therefore skipped entirely.
    """
    if trivy_class == "os-pkgs":
        return "os-package"
    if trivy_type == "jar":
        return "bundled-jar"
    if trivy_class == "lang-pkgs":
        return "conda-requirements-pin"
    return "unknown"


def collect_targets(trivy_path, cve_ids, errors):
    """Resolve each requested CVE ID to its package + fix version(s) from the Trivy JSON."""
    try:
        with open(trivy_path, encoding="utf-8") as handle:
            report = json.load(handle)
    except (OSError, ValueError) as exc:
        errors.append("could not read %s: %s" % (trivy_path, exc))
        return {}

    wanted = {c.upper() for c in cve_ids}
    targets = {}
    for result in report.get("Results") or []:
        for vuln in result.get("Vulnerabilities") or []:
            cve = (vuln.get("VulnerabilityID") or "").upper()
            if cve not in wanted:
                continue
            pkg = vuln.get("PkgName")
            if not pkg:
                continue
            entry = targets.setdefault(
                pkg,
                {
                    "trivy_pkg_name": pkg,
                    "cve_ids": [],
                    "installed_version": vuln.get("InstalledVersion"),
                    "fixed_versions": [],
                    "pkg_paths": [],
                    "trivy_class": result.get("Class"),
                    "trivy_type": result.get("Type"),
                    "fix_shape": classify_fix_shape(result.get("Class"), result.get("Type")),
                    "severities": [],
                },
            )
            if cve not in entry["cve_ids"]:
                entry["cve_ids"].append(cve)
            if vuln.get("Severity") and vuln["Severity"] not in entry["severities"]:
                entry["severities"].append(vuln["Severity"])
            if vuln.get("PkgPath") and vuln["PkgPath"] not in entry["pkg_paths"]:
                entry["pkg_paths"].append(vuln["PkgPath"])
            # Trivy may list several fix branches, e.g. "2.18.8, 2.21.4, 3.1.4".
            for candidate in (vuln.get("FixedVersion") or "").split(","):
                candidate = candidate.strip()
                if candidate and candidate not in entry["fixed_versions"]:
                    entry["fixed_versions"].append(candidate)

    missing = wanted - {c for e in targets.values() for c in e["cve_ids"]}
    for cve in sorted(missing):
        errors.append(
            "%s not found in %s (test mode, or the scan target has diverged)"
            % (cve, trivy_path)
        )
    return targets


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------


def summarize(entry):
    lines = []
    name = entry["trivy_pkg_name"]
    cves = ", ".join(entry["cve_ids"])
    threshold = entry.get("fix_threshold")

    if entry["fix_shape"] == "os-package":
        return [
            "%s (%s): Ubuntu OS package (%s -> %s). Fixed by a base-image / apt update, "
            "NOT by a docker/requirements pin -- conda channels were not queried."
            % (name, cves, entry.get("installed_version"), ", ".join(entry["fixed_versions"]))
        ]

    hit = None
    for channel, data in entry["channels"].items():
        if data.get("arch_ok") and data.get("recommended_floor"):
            hit = (channel, data["recommended_floor"])
            break

    if hit:
        lines.append(
            "%s (%s): %s %s is available on %s for linux-64 + linux-aarch64 -> "
            "floor bump VERIFIED as installable."
            % (name, cves, name, hit[1], hit[0])
        )
    elif any(d.get("name_lookup") != "not_found" for d in entry["channels"].values()):
        lines.append(
            "%s (%s): no build >= %s satisfies both linux-64 and linux-aarch64 -> "
            "treat availability as UNVERIFIED."
            % (name, cves, threshold)
        )
    elif entry["fix_shape"] == "bundled-jar":
        lines.append(
            "%s (%s): a dependency bundled inside a fat JAR, not a conda package. Fixed "
            "only when the upstream tool ships a newer JAR -- see the Java fat-JAR "
            "section of the container-vulns skill." % (name, cves)
        )
    else:
        lines.append(
            "%s (%s): not packaged on conda-forge or bioconda under this name -> a "
            "docker/requirements floor bump does not apply. Check whether it is a "
            "vendored copy or bundled binary instead." % (name, cves)
        )

    for con in entry.get("constrained_by") or []:
        if con["floor_needs_raise"]:
            lines.append(
                "  %s pins '%s'; the current pin '%s' (%s:%s) also admits versions that "
                "exclude %s %s. Raise it to >=%s alongside the bump."
                % (
                    con["package"],
                    con["constraint_at_latest"],
                    con["repo_spec"],
                    con["repo_file"],
                    con["repo_line"],
                    name,
                    threshold,
                    con["min_version_allowing_target"],
                )
            )
        elif con["latest_allows_target"] is False:
            lines.append(
                "  %s pins '%s' and NO available version admits %s %s -- a bump would be "
                "unsolvable until upstream widens the cap."
                % (con["package"], con["constraint_at_latest"], name, threshold)
            )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trivy-json", required=True)
    parser.add_argument(
        "--cve-ids", default="", help="space- and/or comma-separated CVE IDs"
    )
    parser.add_argument("--requirements-dir", default="docker/requirements")
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--no-reverse-scan",
        action="store_true",
        help="skip the reverse dependency-constraint scan",
    )
    args = parser.parse_args()

    cve_ids = [c for c in re.split(r"[\s,]+", args.cve_ids) if c]
    errors = []
    client = ChannelClient(errors)

    targets = collect_targets(args.trivy_json, cve_ids, errors) if cve_ids else {}
    pins = read_requirements(args.requirements_dir)

    packages = []
    for pkg_name in sorted(targets):
        entry = targets[pkg_name]
        fixes = sorted_versions(entry["fixed_versions"]) if entry["fixed_versions"] else []
        # Lowest fix branch at or above what's installed is the floor we actually want.
        installed = entry.get("installed_version")
        usable_fixes = [
            f for f in fixes if not installed or version_ge(f, installed)
        ] or fixes
        entry["fix_threshold"] = usable_fixes[0] if usable_fixes else None

        pin = pins.get(pkg_name.lower())
        entry["repo_requirement"] = pin

        entry["channels"] = {}
        if entry["fix_shape"] != "os-package":
            for channel in CHANNELS:
                entry["channels"][channel] = probe_channel(
                    client, channel, pkg_name, entry["fix_threshold"]
                )

        entry["constrained_by"] = []
        entry["reverse_scan_scope"] = None
        if (
            not args.no_reverse_scan
            and entry["fix_threshold"]
            and entry["fix_shape"] != "os-package"
        ):
            candidates, scoped_file = reverse_scan_candidates(pins, pkg_name)
            entry["reverse_scan_scope"] = {
                "requirements_file": scoped_file,
                "packages_considered": len(candidates),
            }
            entry["constrained_by"] = find_reverse_constraints(
                client, candidates, pkg_name, entry["fix_threshold"], errors
            )

        packages.append(entry)

    summary = []
    for entry in packages:
        summary.extend(summarize(entry))
    if not packages:
        summary.append("No packages resolved from the requested CVE IDs.")

    report = {
        "trivy_json": args.trivy_json,
        "requested_cve_ids": cve_ids,
        "channels_queried": list(CHANNELS),
        "required_subdirs": list(REQUIRED_SUBDIRS),
        "python_build_tag": PYTHON_TAG,
        "summary": summary,
        "packages": packages,
        "probe_errors": errors,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=False)
        handle.write("\n")

    print("Wrote %s" % args.out)
    for line in summary:
        print("  " + line)
    for err in errors:
        print("  probe_error: %s" % err, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
