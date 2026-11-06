# Container Vulnerability Management

Guidance for scanning, triaging, and mitigating container image vulnerabilities
in the viral-ngs Docker image hierarchy.

## Scanning

Container images are scanned for vulnerabilities using [Trivy](https://aquasecurity.github.io/trivy/):

### PR/Push Scans (`docker.yml`)

- Runs on every PR/push after image build
- **SARIF scan**: CRITICAL/HIGH, fixable-only → uploads to GitHub Security tab
- **JSON scan**: CRITICAL/HIGH, fixable-only, **exit-code: 1** → gates CI (fails if CVEs found)
- Both apply `.trivyignore.yaml` and `.trivy-ignore-policy.rego`

### Scheduled Scans (`container-scan.yml`)

- Runs daily at 09:00 UTC on published `main-mega-amd64` image
- **SARIF scan**: CRITICAL/HIGH, fixable-only → uploads to GitHub Security tab
- **JSON scan**: CRITICAL/HIGH/MEDIUM, all vulnerabilities (fixable + non-fixable), **exit-code: 0** → uploaded as artifact for comprehensive debugging and Claude analysis
- **Rebuild first**: when new fixable HIGH/CRITICAL CVEs are detected (i.e. IDs that no open or closed GH issue title names, directly or by a CVE/GHSA alias — `.github/scripts/cve_pipeline.py new-ids`), the `request-rebuild` job dispatches a no-cache rebuild of main's images (`docker.yml` with `force_rebuild=true`, `rescan=true`) instead of triaging. The layer cache otherwise keeps the baseimage's `apt-get upgrade` layer until the upstream micromamba image changes, so Ubuntu security updates published in between never land; a no-cache rebuild picks them up, re-solves every conda requirement, and republishes `main-*` and quay. When it finishes, docker.yml's `request-rescan` job dispatches `container-scan.yml` again with `post_rebuild_run_id` set, which scans the rebuilt image by digest. Most new findings are cleared there and never reach Claude. `cve_pipeline.py route` makes the call; a rescan never routes to another rebuild, so the chain can't loop
- **Triage**: whatever survives the rebuild goes to Claude Sonnet 5.5 on Vertex AI, which triages each one; the workflow files a GitHub issue per CVE (labels: `security`, `cve`). The triage prompt is told the rebuild state (`rebuilt`, `rebuild-failed` or `not-rebuilt`), so "rebuild the image" is never its answer after a successful rebuild. Test mode, `dry_run` and `skip_rebuild` go straight to triage
- **`scan-red` issue, one per rebuild cycle**: the rebuild and its rescan are started by github-actions[bot], and GitHub emails a run's result only to whoever started it, so the `scan-red` job tracks each cycle in an issue (labels: `scan-red`, `security`). It is opened when a scan dispatches a rebuild, and closed by the rescan once the rescan's own findings show the rebuild cleared everything, or every survivor has an issue (linked in the closing comment). It stays open, with a comment, if the rebuild failed, if its `deploy-to-quay` did not succeed (the scan looks at GHCR; users pull from quay.io), if a survivor is still unfiled, or if the rescan's own scan failed (an incomplete scan never closes anything; `gh run rerun <run>` retries it with the same inputs); a later daily scan with nothing new closes a leftover one. Only one CVE rebuild runs at a time: `request-rebuild` skips the dispatch while another CVE-pipeline dispatch of docker.yml (named by its `run-name`) is in flight (a concurrency group keeps two scans from both dispatching). Decisions live in `cve_pipeline.py scan-red`
- **Open `cve` issues are re-checked every run** (`cve_pipeline.py backlog`, acted on by the `update-cve-issues` job). Only issues this workflow filed count (bot-authored, not labeled `test`). An issue whose advisory, by ID or CVE/GHSA alias, is no longer in the scan is closed with a comment linking the run, and reopened if the advisory comes back (the scan never files an advisory twice), unless someone closed it again since: a fix PR's or a maintainer's close stands. A report with no Ubuntu result closes nothing, and nothing is closed or reopened once `main-mega-amd64` no longer points at the digest the scan saw (a push landed in between; the next daily scan decides). An issue whose findings are all Ubuntu packages is waiting on Ubuntu: the fix-PR agent never sees it, and each run asks apt, in the image `Dockerfile.baseimage` builds FROM, whether the archive's candidate version has reached Trivy's fixed version. Once it has, the route is `rebuild` for that issue too. The scan passes the `<package>=<version>` pairs it saw published through the rebuild (`pre_rebuild_versions`). If a successful rebuild's rescan still finds the advisory, each of those versions gets a hidden `cve-retry:<package>=<version>` marker on the issue and a comment asking for a maintainer (apt may be holding the package back), and is never retried again; a newer version, including one Ubuntu published while the rebuild ran, gets its own rebuild
- **Auto-fix PR (`cve-fix-pr` job)**: after new issues are filed, Claude assesses *all* currently-open `cve` issues, except the Ubuntu-only and no-longer-reported ones above, and either opens a single `cve-fix` + `security` PR — only if every open issue meets the moderate confidence bar (version-floor bump, precedent-mirroring Dockerfile mitigation, or precedent-mirroring `.trivyignore.yaml` entry; all-or-nothing) — or uploads a `SKIPPED.md` decision-log artifact. The PR is authored with a GitHub App token so it triggers the required CI checks. Its title, which is what release notes list, names the advisories it closes: the agent writes `TITLE.txt`, and `cve_pipeline.py pr-title` uses it only if it is plain text naming every closed advisory, otherwise it lists the closed issues' IDs grouped by package.
- **Recovery / manual re-run**: the auto-fix job only fires when the scan finds a *new* CVE, so once an issue is filed the CVE is no longer "new" and a failed or skipped fix run will not retry on its own. To re-run it on demand against all currently-open `cve` issues, dispatch manually: `gh workflow run container-scan.yml -f force_fix_pr=true`. To rebuild main's images without cache by hand: `gh workflow run docker.yml --ref main -f force_rebuild=true` (add `-f rescan=true` to have the CVE pipeline rescan and triage afterwards; docker.yml refuses a dispatch on any ref but main). To send new CVEs straight to Claude without a rebuild: `gh workflow run container-scan.yml -f skip_rebuild=true`
- **Conda availability probe**: both Claude steps are preceded by `.github/scripts/probe-conda-availability.py`, which queries the conda-forge and bioconda APIs and writes `.cve-fix-context/conda-availability.json`. Per affected package it reports whether the Trivy `FixedVersion` exists on a channel, whether it is built for **both** `linux-64` and `linux-aarch64` (or is `noarch`), whether a requirements pin can fix the finding at all (`fix_shape` — an Ubuntu `os-package` or a `bundled-jar` cannot), and which other pinned packages carry a constraint that would exclude the fix version (`constrained_by`, e.g. `pyopenssl<50` capping `cryptography`). Uploaded as a workflow artifact. Before this existed the agents had no network tool at all and skipped *every* version-floor bump as unverifiable — the primary category the pipeline was built for.
- **CI verifies ARM64, so the agent does not have to**: `docker.yml` builds all six tiers for both amd64 and native arm64 on every PR and re-scans each with Trivy, so a solver conflict from a floor bump fails required checks before merge. The auto-fix prompt says so explicitly and instructs the agent to reserve SKIP for "I cannot determine *what* the fix is", not "I cannot prove the fix works". Unproven solvability used to be the standard reason for skipping.
- **Steering the bot via issue comments**: maintainer comments on a `cve` issue are passed to the auto-fix agent and outrank the machine-generated issue body — use them to specify or amend a fix. Comments are hard-filtered to `OWNER`/`MEMBER`/`COLLABORATOR` in the workflow's `jq`, not in the prompt: this is a public repo, so an unfiltered comments array would be an untrusted-instruction channel into an agent that edits files and opens PRs.
- **Agent guardrails**: each agent runs in a job of its own (`triage`, `cve-fix-pr`) whose `GITHUB_TOKEN` is read-only. Writes happen in jobs no agent ran in: `scan` uploads SARIF, `file-issues` files the triage reports (only for the IDs it asked about), and `cve-fix-pr-open` applies the auto-fix agent's patch to a fresh checkout and is the only job holding the GitHub App token. That split matters because claude-code-action exports its `github_token` into the agent's environment, and an agent with Write has to be assumed able to run code in its job (`git grep -O<cmd>` is on its allowlist; `.git/config` it writes is obeyed by any later git command). `cve-fix-pr-open` treats the patch as untrusted: it must stay inside the path allowlist (`docker/requirements/*.txt`, `docker/Dockerfile.*`, `.trivyignore.yaml`), change only the text content of regular files (no binary content, symlinks, modes or renames), and pass gitleaks along with the PR body and title. It also refuses any Dockerfile hunk that changes a line within 3 lines of an `apt-get` command (see "Never append package names to the apt line" below). It checks out `github.sha` itself rather than taking a commit from the agent's job. The auto-fix agent has domain-scoped `WebFetch` (no `WebSearch`, since its output is committed); the triage agent gets `WebSearch` too.
- The Vertex/WIF infra used here is documented in `.agents/skills/claude-on-vertex-ci/SKILL.md`

**Key difference**: The scheduled scan's JSON output includes MEDIUM severity and non-fixable vulnerabilities to provide full visibility for risk assessment and debugging, while the control flow (issue filing, CI gates) still filters to fixable HIGH/CRITICAL only.

Per-CVE exceptions go in `.trivyignore.yaml` with mandatory justification comments.

The file uses Trivy's **YAML** ignore format rather than the one-CVE-per-line plain text
format, because a plain-text line matches a CVE ID *globally* -- it suppresses that CVE
anywhere in the image, including copies nobody assessed. Every entry must be scoped:

- `paths: [...]` -- match by file location; globs work (`opt/conda/share/snpeff-*/snpEff.jar`),
  which survives version bumps. Use for findings that carry a `PkgPath`, e.g. a dependency
  inside a fat JAR.
- `purls: [...]` -- match by package identity + version (`pkg:pypi/setuptools@70.3.0`). Use
  for findings with **no** `PkgPath`. In particular, anything Trivy derives from an embedded
  SBOM (the record shows `"AnalyzedBy": "sbom"`) has no file path, so `paths` cannot match it
  -- pip ships such an SBOM at `pip/_vendor/bom.cdx.json` describing its vendored deps.

Scoping means that if the same CVE later appears somewhere unassessed, it is reported
instead of silently suppressed.

**Key entries by CVE ID when there is one.** Trivy reports a GitHub advisory under its GHSA ID
only while no CVE is assigned. Once one is, the finding's `VulnerabilityID` becomes the CVE and
the GHSA ID moves to `VendorIDs`. Ignore entries match `VulnerabilityID` only, so a GHSA-keyed
entry stops matching and the finding comes back. When the scheduled scan sees a finding whose
alias is already tracked by an issue, it files no new issue and logs a warning naming both IDs:
that is the cue to re-key the entry.

## Rego Policy (`.trivy-ignore-policy.rego`)

The Rego policy filters CVEs that are architecturally inapplicable to ephemeral batch containers:

- **AV:P** (Physical access required) -- containers are cloud-hosted
- **AV:A** (Adjacent network required) -- no attacker on same network segment
- **AV:L + UI:R** (Local + user interaction) -- no interactive sessions
- **AV:L + PR:H** (Local + high privileges) -- containers run non-root
- **AV:L + S:U** (Local + scope unchanged) -- attacker already has code execution and impact stays within the ephemeral container

Changes to this policy should be reviewed carefully. The comments in the file explain the rationale and risk for each rule.

## Common Vulnerability Sources

**Python transitive deps**: Pin minimum versions in `docker/requirements/*.txt`. Prefer conda packages over pip. Check conda-forge availability before assuming a version exists -- conda-forge often lags PyPI by days/weeks.

**Java fat JARs** (picard, gatk, snpeff, fgbio): Bioinformatics Java tools are distributed as uber JARs with all dependencies bundled inside. Trivy detects vulnerable libraries (log4j, commons-compress, etc.) baked into these JARs. Version bumps can cause ARM64 conda solver conflicts because Java tools pull in openjdk -> harfbuzz -> icu version chains that clash with other packages (r-base, boost-cpp, pyicu). Always check:
1. Whether the tool is actually flagged by Trivy (don't bump versions unnecessarily)
2. Whether the CVE applies (e.g., log4j 1.x is NOT vulnerable to Log4Shell)
3. Whether the desired version resolves on ARM64 before pushing

**Go binaries**: Some conda packages bundle compiled Go binaries (e.g., mafft's `dash_client`, google-cloud-sdk's `gcloud-crc32c`). If the binary is unused, delete it in the Dockerfile. Delete from **both** the installed location and `/opt/conda/pkgs/*/` (conda package cache) -- Trivy scans the full filesystem.

**Vendored copies**: Packages like google-cloud-sdk and setuptools bundle their own copies of Python libraries that may be older than what's in the conda environment. Trivy flags these vendored copies separately. Options: delete the vendored directory (if not needed at runtime), or accept the risk in `.trivyignore.yaml` with justification. Prefer an ignore entry over deleting an embedded SBOM: the SBOM usually describes several vendored packages, so deleting it to silence one mis-mapped CVE hides real future findings in the rest.

## ARM64 Solver Conflicts

The conda solver on ARM64 (linux-aarch64) is more constrained than amd64 because fewer package builds exist. Common conflict patterns:

- **icu version conflicts**: Many packages (openjdk, r-base, boost-cpp, pyicu) pin specific icu version ranges. Bumping one package can make the entire environment unsolvable.
- **libdeflate/htslib conflicts**: lofreq 2.1.5 pins old htslib/libdeflate versions that conflict with newer pillow/libtiff.
- **openjdk version escalation**: snpeff 5.2+ requires openjdk>=11, 5.3+ requires openjdk>=21. Higher openjdk versions pull in harfbuzz->icu chains that conflict with everything.

When a solver conflict occurs: revert the change, check what version the solver was picking before, and pin to that exact version if it already addresses the CVE.

## Never append package names to the apt line

`docker/Dockerfile.baseimage` runs a plain `apt-get upgrade -y`. Until November 2026 the CVE
bot "fixed" Ubuntu-package CVEs by appending package names to that line (gpgv, openssl, curl,
libbrotli1, apt: commits 423d1625, fbab14c4, a1dfe313, 2823b591, 3cf00139). The names did
nothing: `apt-get upgrade <names>` upgrades every installed package anyway. Each fix worked
only because changing the line's text invalidated the cached layer. A no-cache rebuild
(`force_rebuild`) does the same thing honestly, with no commit, and the pipeline now runs one
before triage. An `os-package` finding that survives that rebuild is waiting on Ubuntu to
publish the fix. Its issue stays open, out of the fix-PR agent's way, and the daily scan
dispatches another rebuild once apt sees the fixed version in the archive.

## Mitigation Decision Process

When triaging a CVE:

0. **Does a no-cache rebuild clear it?** The pipeline dispatches one before triage. If the
   finding survived (or you are triaging by hand), an Ubuntu `os-package` finding usually
   means the fixed version is not in the archive yet, and a conda finding means no existing
   floor allows the fix.
1. **Check the CVSS vector** -- does the Rego policy already filter it?
2. **Identify the source package** -- use Trivy JSON output (`PkgName`, `PkgPath`, `InstalledVersion`)
3. **Check if a fix version exists on conda-forge/bioconda** -- not just on PyPI. Run
   `.github/scripts/probe-conda-availability.py --trivy-json trivy-results.json --cve-ids "<CVE-ID>" --out /tmp/avail.json`
   (works locally against a downloaded scan artifact) rather than checking by hand; it also
   reports which other pinned packages constrain the version you want.
4. **Test on ARM64** -- solver conflicts are the most common failure mode. Opening a PR is
   a legitimate way to test this: `docker.yml` builds native arm64 for every tier.
5. **If the fix version conflicts**: consider whether the CVE is exploitable in your deployment model. Document the risk assessment in `.trivyignore.yaml` or `vulnerability-mitigation-status.md`.
6. **If the vulnerable code is unused**: delete the binary/file inline in the Dockerfile (same RUN layer as install to avoid bloating images)

## Key Files

| File | Purpose |
|------|---------|
| `.trivy-ignore-policy.rego` | Rego policy for class-level CVE filtering |
| `.trivyignore.yaml` | Per-CVE exceptions, scoped by `paths` or `purls`, with justifications |
| `.github/workflows/docker.yml` | Build-time scanning (SARIF + JSON), amd64 + native arm64; `workflow_dispatch` with `force_rebuild` for a no-cache rebuild of main's images |
| `.github/workflows/container-scan.yml` | Daily scheduled scanning, rebuild-first routing, `scan-red` issue, Claude triage, auto-fix PR |
| `.github/scripts/cve_pipeline.py` | The pipeline's decisions: `new-ids`, `route`, `scan-red`, `backlog`, `backlog-actions`, `pr-title` (unit-tested in `test_cve_pipeline.py`) |
| `.github/scripts/probe-conda-availability.py` | Pre-agent conda-forge/bioconda availability + reverse-constraint probe |
| `vulnerability-mitigation-status.md` | Local-only tracking doc (not committed) |
