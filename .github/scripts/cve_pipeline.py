#!/usr/bin/env python3
"""Advisory-ID decisions for the container-scan.yml CVE pipeline.

The workflow fetches the inputs (Trivy report, GitHub issues); this script only decides,
so the matching rules are unit-tested instead of living in inline bash.

  new-ids   Fixable HIGH/CRITICAL findings in a Trivy report that no issue tracks yet.
            An issue tracks a finding when its title names the finding's VulnerabilityID
            or one of its CVE/GHSA aliases in VendorIDs. Trivy reports a GitHub advisory
            under its CVE once one is assigned and moves the GHSA ID to VendorIDs, so
            without the alias check a GHSA-titled issue would be filed again under the CVE.
  names     Exit 0 if --text names --id as a whole ID, 1 if not: the file-issues job's
            check that an issue title will be recognized by new-ids.
  pr-title  Title for the auto-fix PR. Uses the agent's TITLE.txt when it is one line of
            plain text, at most MAX_TITLE_LENGTH characters, and names every advisory the
            PR closes; otherwise builds one from the closed issues' titles.

IDs match as whole IDs only: CVE-2026-1123 is not named by a title about CVE-2026-11234.

Usage:
    cve_pipeline.py new-ids --trivy-json trivy-results.json \\
        --issues-json all-issues.json --out new-ids.txt
    cve_pipeline.py names --id CVE-2026-1123 --text "$title"
    cve_pipeline.py pr-title --title-file .claude-fix-pr/TITLE.txt \\
        --pr-body .claude-fix-pr/PROCEED.md \\
        --issues-json .cve-fix-context/open-issues.json --out pr-title.txt

new-ids and pr-title write their result to --out and report decisions on stdout, as
workflow annotations when run under GitHub Actions.
"""

import argparse
import json
import os
import re
import string
import sys

ADVISORY_ID = re.compile(r'CVE-\d{4}-\d{4,}|GHSA(?:-[0-9a-z]{4}){3}', re.IGNORECASE)
# GitHub's closing keywords ("Closes #12", "fixes #12", "Resolved: #12").
CLOSING_REF = re.compile(r'\b(?:close[sd]?|fix(?:e[sd])?|resolve[sd]?):?\s+#(\d+)\b', re.IGNORECASE)
# The triage job titles issues "[<ID>] <package>: <description>". Maven packages contain
# ':' themselves, so the package ends at the first ': '.
ISSUE_TITLE = re.compile(r'^\[[^\]]+\]\s*(.+?):\s')

MAX_TITLE_LENGTH = 200  # GitHub allows 256
# Plain text plus what version specifiers need. No '@' (mentions), '#' (issue links),
# backticks, brackets or other markdown: PR titles are rendered in release notes.
TITLE_CHARS = frozenset(string.ascii_letters + string.digits + " .,:;()<>=+/_'~-")
GENERIC_TITLE = 'Auto-fix: address open CVE issues'


def annotate(out, level, message):
    if os.environ.get('GITHUB_ACTIONS') == 'true':
        print(f'::{level}::{message}', file=out)
    else:
        print(f'{level.upper()}: {message}', file=out)


def mentions(text, advisory_id):
    """True if text names advisory_id as a whole ID."""
    pattern = rf'(?<![0-9A-Za-z]){re.escape(advisory_id)}(?![0-9A-Za-z])'
    return re.search(pattern, text, re.IGNORECASE) is not None


def tracking_issue(advisory_id, issues):
    return next((i for i in issues if mentions(i['title'], advisory_id)), None)


def fixable_findings(report):
    for result in report.get('Results') or []:
        for vuln in result.get('Vulnerabilities') or []:
            if vuln.get('Severity') in ('HIGH', 'CRITICAL') and vuln.get('FixedVersion'):
                yield vuln


def new_advisories(report, issues, out=sys.stdout):
    """Sorted IDs of fixable HIGH/CRITICAL findings that no issue title names, by ID or alias."""
    aliases = {}
    for vuln in fixable_findings(report):
        vuln_id = vuln['VulnerabilityID']
        aliases.setdefault(vuln_id, set()).update(
            a for a in vuln.get('VendorIDs') or [] if ADVISORY_ID.fullmatch(a) and a != vuln_id)
    if not aliases:
        annotate(out, 'notice', 'No fixable HIGH/CRITICAL CVEs in scan.')
        return []

    new = []
    for vuln_id in sorted(aliases):
        tracked = tracking_issue(vuln_id, issues)
        if tracked:
            print(f'  existing issue #{tracked["number"]} for {vuln_id}, skipping', file=out)
            continue
        for alias in sorted(aliases[vuln_id]):
            tracked = tracking_issue(alias, issues)
            if tracked:
                annotate(out, 'warning',
                         f'{vuln_id} is already tracked as {alias} in #{tracked["number"]}, so no new issue '
                         f'is filed. Trivy now reports it as {vuln_id} and matches .trivyignore.yaml entries '
                         f'on that ID only: any entry keyed {alias} must be re-keyed to {vuln_id}.')
                break
        else:
            print(f'  NEW: {vuln_id}', file=out)
            new.append(vuln_id)
    annotate(out, 'notice', f'Found {len(new)} new fixable CVE(s)')
    return new


def issue_package(title):
    """The package from a "[<ID>] <package>: ..." title, or None if absent or not plain text."""
    m = ISSUE_TITLE.match(title)
    if m and set(m.group(1)) <= TITLE_CHARS and '://' not in m.group(1):
        return m.group(1).strip()
    return None


def closed_advisories(pr_body, issues, out):
    """{advisory ID: package} for each open cve issue the PR body closes, in closing order."""
    by_number = {i['number']: i for i in issues}
    advisories = {}
    for number in dict.fromkeys(int(n) for n in CLOSING_REF.findall(pr_body)):
        closed = by_number.get(number)
        if closed is None:
            annotate(out, 'warning', f'The PR body closes #{number}, which is not an open cve issue')
            continue
        m = ADVISORY_ID.search(closed['title'])
        if not m:
            annotate(out, 'warning', f'#{number} has no advisory ID in its title: {closed["title"]!r}')
            continue
        advisories.setdefault(m.group(0), issue_package(closed['title']))
    return advisories


def rejection_reason(agent_title, advisory_ids):
    """Why agent_title can't be the PR title, or None if it can."""
    if agent_title is None:
        return 'the agent wrote no TITLE.txt'
    lines = [line.strip() for line in agent_title.splitlines() if line.strip()]
    if not lines:
        return 'TITLE.txt is empty'
    if len(lines) > 1:
        return f'TITLE.txt has {len(lines)} lines'
    title = lines[0]
    if len(title) > MAX_TITLE_LENGTH:
        return f'it is {len(title)} characters (limit {MAX_TITLE_LENGTH})'
    disallowed = sorted(set(title) - TITLE_CHARS)
    if disallowed:
        return f'it contains disallowed characters {"".join(disallowed)!r}'
    if '://' in title:
        return 'it contains a URL'
    missing = [a for a in advisory_ids if not mentions(title, a)]
    if missing:
        return f'it does not name {", ".join(missing)}'
    return None


def fallback_title(advisories):
    """A title built from the closed issues: IDs grouped by package, shortened to fit.

    Package names go first when shortening: the IDs are what release notes and
    `git log --grep CVE-` need.
    """
    groups = {}
    for advisory_id, package in advisories.items():
        groups.setdefault(package, []).append(advisory_id)
    parts = [' '.join(filter(None, (package, ', '.join(ids)))) for package, ids in groups.items()]
    title = 'Auto-fix: ' + '; '.join(parts)
    if len(title) <= MAX_TITLE_LENGTH:
        return title
    title = 'Auto-fix: ' + ', '.join(advisories)
    if len(title) <= MAX_TITLE_LENGTH:
        return title
    packages = [p for p in groups if p]
    count = f'{len(advisories)} advisories'
    title = f'Auto-fix: {", ".join(packages)} ({count})' if packages else f'Auto-fix: {count}'
    if len(title) <= MAX_TITLE_LENGTH:
        return title
    return f'Auto-fix: {count} in {len(packages)} packages'


def pr_title(agent_title, pr_body, issues, out=sys.stdout):
    advisories = closed_advisories(pr_body, issues, out)
    if not advisories:
        annotate(out, 'warning', 'The PR body closes no open cve issue; using a generic PR title')
        return GENERIC_TITLE
    reason = rejection_reason(agent_title, list(advisories))
    if reason is None:
        title = agent_title.strip()
        print(f"Using the agent's PR title: {title}", file=out)
        return title
    title = fallback_title(advisories)
    annotate(out, 'warning', f"Not using the agent's PR title ({reason}); using {title!r} instead")
    return title


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('new-ids', help='fixable findings no issue tracks yet')
    p.add_argument('--trivy-json', required=True)
    p.add_argument('--issues-json', required=True, help='every issue, open and closed: [{number, title}]')
    p.add_argument('--out', required=True, help='written with the new IDs, space-separated')
    p = commands.add_parser('names', help='exit 0 if --text names --id as a whole ID')
    p.add_argument('--id', required=True, dest='advisory_id')
    p.add_argument('--text', required=True)
    p = commands.add_parser('pr-title', help='title for the auto-fix PR')
    p.add_argument('--title-file', required=True, help="the agent's TITLE.txt (may be absent)")
    p.add_argument('--pr-body', required=True)
    p.add_argument('--issues-json', required=True, help='open cve issues: [{number, title, ...}]')
    p.add_argument('--out', required=True, help='written with the title')
    args = parser.parse_args(argv)

    if args.command == 'names':
        return 0 if mentions(args.text, args.advisory_id) else 1
    with open(args.issues_json) as f:
        issues = json.load(f)
    if args.command == 'new-ids':
        with open(args.trivy_json) as f:
            result = ' '.join(new_advisories(json.load(f), issues, sys.stdout))
    else:
        try:
            with open(args.title_file) as f:
                agent_title = f.read()
        except FileNotFoundError:
            agent_title = None
        with open(args.pr_body) as f:
            result = pr_title(agent_title, f.read(), issues, sys.stdout)
    with open(args.out, 'w') as f:
        f.write(result + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
