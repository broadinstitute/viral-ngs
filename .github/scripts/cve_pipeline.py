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
  route     Where a scan's new IDs go: `rebuild` (dispatch a no-cache rebuild of the main
            images first; its rescan triages whatever survives), `triage` (straight to
            Claude) or `none`. A rebuild is also due when Ubuntu has published the fix an
            open os-package issue is waiting for (backlog's retry IDs). A post-rebuild run
            never routes to `rebuild`, which is what keeps the scan -> rebuild -> rescan
            chain from looping.
  scan-red  What to do with the `scan-red` issue that tracks one rebuild cycle: open it
            when the daily scan dispatches a rebuild, close it when the rescan's findings
            show the rebuild cleared everything or every survivor has an issue, and leave it
            open, with a comment, when something needs a maintainer: the rebuild failed or
            never reached quay.io, a survivor is unfiled, or the rescan itself failed
            (--scan-result; an incomplete scan never closes anything). Writes
            {action: open|comment|close|none, number, title, body, summary} as JSON.
  backlog   What this run's Trivy report says about the open `cve` issues this workflow
            filed (bot-authored, not labeled `test`): `gone` when no finding has the
            title's ID or an alias of it, so the issue can be closed (only from a report
            that scanned the OS); `waiting` when every finding is an Ubuntu package, which
            no PR fixes, so the fix-PR agent never sees it; `retry` for a waiting issue
            once the apt probe (--published) shows Ubuntu's archive has a fixed version
            that no completed rebuild has already failed to install; `reopen` for a closed
            issue (--closed-issues-json) whose latest close was the scan's own close as gone,
            once its advisory is back. Also records the digest the report scanned and the
            retries' "<package>=<version>" pairs, which request-rebuild passes through the
            rebuild to its rescan. Writes JSON, and with --waiting-pkgs-out the
            "<package> <fixed version>" lines the probe checks.
  backlog-actions
            The issue updates backlog calls for, as JSON [{action: close|reopen|comment,
            number, body}]: close each gone issue and reopen each returning one, but only
            while main-mega-amd64 (--main-digest) is still the image the scan saw; and after
            a successful rebuild's rescan, mark each published version that rebuild was
            dispatched for (--pre-rebuild-versions) but did not install as tried (a
            `cve-retry` marker) and ask for a maintainer.

IDs match as whole IDs only: CVE-2026-1123 is not named by a title about CVE-2026-11234.

Usage:
    cve_pipeline.py new-ids --trivy-json trivy-results.json \\
        --issues-json all-issues.json --out new-ids.txt
    cve_pipeline.py names --id CVE-2026-1123 --text "$title"
    cve_pipeline.py pr-title --title-file .claude-fix-pr/TITLE.txt \\
        --pr-body .claude-fix-pr/PROCEED.md \\
        --issues-json .cve-fix-context/open-issues.json --out pr-title.txt
    cve_pipeline.py route --cve-ids "$CVE_IDS" --test-mode false --dry-run false \\
        --skip-rebuild false --post-rebuild-run-id "$RUN_ID" --out route.txt
    cve_pipeline.py scan-red --route rebuild --cve-ids "$CVE_IDS" --dispatch dispatched \\
        --open-issue-json open-scan-red.json --issues-json all-issues.json \\
        --trivy-json trivy-results.json --repo-url "$REPO_URL" --run-id "$GITHUB_RUN_ID" \\
        --image "$IMAGE_REF" --date 2026-11-05 --out decision.json
    cve_pipeline.py backlog --trivy-json trivy-results.json --issues-json open-cve-issues.json \
        [--closed-issues-json closed-cve-issues.json] [--published published.txt] \
        [--waiting-pkgs-out waiting-pkgs.txt] --out backlog.json
    cve_pipeline.py backlog-actions --backlog backlog.json --rebuild-state rebuilt \
        --post-rebuild-run-id "$RUN_ID" --pre-rebuild-versions "$PRE_REBUILD_VERSIONS" \
        --main-digest "$MAIN_DIGEST" --repo-url "$REPO_URL" --run-id "$GITHUB_RUN_ID" \
        --image "$IMAGE_REF" --out actions.json

new-ids, pr-title, route, scan-red, backlog and backlog-actions write their result to --out and report decisions
on stdout, as workflow annotations when run under GitHub Actions.
"""

import argparse
import datetime
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


def advisory_aliases(report):
    """{VulnerabilityID: its CVE/GHSA aliases} for the fixable HIGH/CRITICAL findings in a report."""
    aliases = {}
    for vuln in fixable_findings(report):
        vuln_id = vuln['VulnerabilityID']
        aliases.setdefault(vuln_id, set()).update(
            a for a in vuln.get('VendorIDs') or [] if ADVISORY_ID.fullmatch(a) and a != vuln_id)
    return aliases


def new_advisories(report, issues, out=sys.stdout):
    """Sorted IDs of fixable HIGH/CRITICAL findings that no issue title names, by ID or alias."""
    aliases = advisory_aliases(report)
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


def route(cve_ids, test_mode=False, dry_run=False, skip_rebuild=False, post_rebuild=False, retry_ids=()):
    """Where a scan's new IDs go: 'rebuild', 'triage' or 'none'.

    A no-cache rebuild runs `apt-get upgrade -y` against today's Ubuntu archive and
    re-solves every conda floor, which clears most new findings without a PR. So new IDs
    go to a rebuild first, except when a rebuild already happened (post_rebuild), when a
    maintainer asked for triage directly (skip_rebuild), and in test or dry runs, which
    must not republish main's images. retry_ids (backlog's) are already filed, so they
    route to a rebuild under the same conditions and otherwise nowhere.
    """
    if test_mode:
        return 'triage'
    no_rebuild = dry_run or skip_rebuild or post_rebuild
    if cve_ids:
        return 'triage' if no_rebuild else 'rebuild'
    if retry_ids and not no_rebuild:
        return 'rebuild'
    return 'none'


SCAN_RED_MARKER = '<!-- scan-red:{} -->'
SCAN_RED_MARKER_RE = re.compile(r'<!-- scan-red:(\S+?) -->')
DISPATCH_TEXT = {
    'dispatched': ('Dispatched a no-cache rebuild of the main images (`docker.yml`, `force_rebuild=true`). '
                   'Its rescan triages whatever the rebuild does not clear.'),
    'in-flight': ('A no-cache rebuild of the main images (`docker.yml`, `force_rebuild=true`) is already in '
                  'flight, so no second one was dispatched. Its rescan covers these too.'),
}
DISPATCH_FAILED_TEXT = ('The no-cache rebuild could not be dispatched (see the `request-rebuild` job of this '
                        'run). The next daily scan tries again.')
LIFECYCLE_TEXT = ('The rebuild\'s rescan closes this issue when the rebuild clears these, or once Claude triage '
                  'has filed the survivors as `cve` issues. It stays open, with a comment, if something needs a '
                  'maintainer.')
RERUN_TEXT = 'To retry: `gh workflow run docker.yml --ref main -f force_rebuild=true -f rescan=true`.'


def _markers(issue):
    """The scan-red markers in an issue's body and comments, oldest first."""
    texts = [issue.get('body') or ''] + [c.get('body') or '' for c in issue.get('comments') or []]
    return [m for text in texts for m in SCAN_RED_MARKER_RE.findall(text)]


def _cell(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def findings_table(report, advisory_ids):
    """A markdown table of the fixable HIGH/CRITICAL findings for advisory_ids."""
    groups = {}
    for result in (report or {}).get('Results') or []:
        kind = result.get('Type') or result.get('Class') or ''
        for vuln in result.get('Vulnerabilities') or []:
            if (vuln.get('VulnerabilityID') not in advisory_ids
                    or vuln.get('Severity') not in ('HIGH', 'CRITICAL') or not vuln.get('FixedVersion')):
                continue
            key = (vuln['VulnerabilityID'], vuln['Severity'], kind,
                   vuln.get('InstalledVersion') or '?', vuln['FixedVersion'])
            packages = groups.setdefault(key, [])
            if vuln.get('PkgName', '?') not in packages:
                packages.append(vuln.get('PkgName', '?'))
    rows = ['| Advisory | Severity | Packages | Type | Installed → fixed |', '|---|---|---|---|---|']
    for (vuln_id, severity, kind, installed, fixed), packages in sorted(groups.items()):
        rows.append(f'| {_cell(vuln_id)} | {severity} | {_cell(", ".join(packages))} | {_cell(kind)} '
                    f'| {_cell(installed)} → {_cell(fixed)} |')
    for vuln_id in sorted(set(advisory_ids) - {key[0] for key in groups}):
        rows.append(f'| {_cell(vuln_id)} | | | | not in the Trivy report |')
    return '\n'.join(rows)


def _advisories(n):
    return f'{n} new fixable HIGH/CRITICAL {"advisory" if n == 1 else "advisories"}'


def scan_red(*, route, cve_ids, retry_ids=(), pre_rebuild_ids=(), post_rebuild_run_id='', rebuild_result='',
             deploy_result='', dispatch='', file_issues_result='', open_issue=None, issues=(), report=None,
             repo_url='', run_id='', image='', date='', scan_result='success'):
    """How to update the scan-red issue that tracks one rebuild cycle. See the module doc."""
    ids = sorted(set(cve_ids))
    number = open_issue['number'] if open_issue else None
    markers = _markers(open_issue) if open_issue else []
    run_link = f'[run {run_id}]({repo_url}/actions/runs/{run_id})'

    def decide(action, body='', title=None):
        return {'action': action, 'number': number if action in ('comment', 'close') else None,
                'title': title, 'body': body, 'summary': body}

    if scan_result != 'success':
        # An incomplete scan says nothing about the image, so it never closes a cycle. A failed
        # daily scan leaves everything to the next one; a failed rescan would otherwise leave
        # the cycle's issue claiming a rebuild is in progress, and nobody gets its run's email.
        if not post_rebuild_run_id:
            return decide('none')
        rebuild_link = f'[run {post_rebuild_run_id}]({repo_url}/actions/runs/{post_rebuild_run_id})'
        marker = f'rescan-failed:{run_id}'
        body = (f'The rescan in {run_link}, after the no-cache rebuild in {rebuild_link}, did not complete '
                f'(`{scan_result or "unknown"}`), so nothing was triaged or filed, and whether the rebuild fixed '
                f'anything is unknown. The rebuild itself: build `{rebuild_result or "unknown"}`, deploy to '
                f'quay.io `{deploy_result or "unknown"}`. The next daily scan checks main\'s image again.\n\n'
                f'To retry the rescan now: `gh run rerun {run_id}`. **If it fails again, this needs a '
                f'maintainer.**\n\n' + SCAN_RED_MARKER.format(marker))
        if open_issue:
            return decide('none', body) if marker in markers else decide('comment', body)
        return decide('open', body, f'Rescan after a no-cache rebuild failed (since {date})')

    if post_rebuild_run_id:
        rebuild_link = f'[run {post_rebuild_run_id}]({repo_url}/actions/runs/{post_rebuild_run_id})'
        if rebuild_result != 'success':
            marker = f'rebuild-failed:{post_rebuild_run_id}'
            body = (f'The no-cache rebuild in {rebuild_link} did not succeed (`{rebuild_result or "unknown"}`), '
                    f'so the main images were not republished. {run_link} scanned `{image}` and found '
                    f'{_advisories(len(ids))}; any it found go to Claude triage, which cannot tell whether a '
                    f'rebuild would have cleared them.\n\n**This needs a maintainer**: fix the build, then '
                    f'retry the rebuild. {RERUN_TEXT}\n\n'
                    + (findings_table(report, ids) + '\n\n' if ids else '')
                    + SCAN_RED_MARKER.format(marker))
            if open_issue:
                return decide('none', body) if marker in markers else decide('comment', body)
            return decide('open', body, f'No-cache rebuild of the main images failed (since {date})')

        # What the rebuild cleared comes from the rescan's findings, not from the new IDs: an
        # advisory that got an issue while the rebuild ran is not new, yet still in the image.
        if report is None:
            aliases, present = {}, None
        else:
            aliases = advisory_aliases(report)
            present = {}  # every ID and alias of a finding -> its VulnerabilityID
            for vuln_id, names in aliases.items():
                present[vuln_id] = vuln_id
                for name in names:
                    present.setdefault(name, vuln_id)
        if present is None:
            cleared, lingering, survivors = [], [], ids
        else:
            # Gone from the report altogether, or still in it but no longer fixable HIGH/CRITICAL.
            anywhere = _advisory_names(report)
            gone = [i for i in sorted(set(pre_rebuild_ids)) if i not in present]
            cleared = [i for i in gone if i.upper() not in anywhere]
            lingering = [i for i in gone if i.upper() in anywhere]
            survivors = sorted(set(ids) | {present[i] for i in pre_rebuild_ids if i in present})

        def tracked_by(vuln_id):
            for name in [vuln_id, *sorted(aliases.get(vuln_id, ()))]:
                tracked = tracking_issue(name, issues)
                if tracked:
                    return tracked
            return None

        tracking = {v: tracked_by(v) for v in survivors}
        untracked = [v for v in survivors if not tracking[v]]

        # Reasons to keep the cycle open: (marker, paragraph, title if this opens the issue).
        problems = []
        if deploy_result != 'success':
            problems.append((
                f'deploy-failed:{post_rebuild_run_id}',
                f"The rebuild's `deploy-to-quay` did not succeed (`{deploy_result or 'unknown'}`), so quay.io "
                f'still serves the previous images, though the GHCR image scanned here was rebuilt. '
                f'**This needs a maintainer**: see {rebuild_link}. {RERUN_TEXT}',
                f'No-cache rebuild did not reach quay.io (since {date})'))
        if untracked:
            problems.append((
                f'not-filed:{run_id}',
                f'{", ".join(untracked)} survived the rebuild but {"has" if len(untracked) == 1 else "have"} '
                f'no issue yet (`file-issues` in {run_link}: `{file_issues_result or "unknown"}`). Still '
                f'unfiled, so the next daily scan picks {"it" if len(untracked) == 1 else "them"} up again.\n\n'
                + findings_table(report, untracked),
                f'Advisories survived the no-cache rebuild but were not filed (since {date})'))
        if present is None:
            problems.append((
                f'no-report:{run_id}',
                f"The rescan's Trivy report was not available to this job, so what the rebuild cleared is "
                f'unknown. The next daily scan closes this issue if nothing new remains.',
                f'Rescan after a no-cache rebuild could not be checked (since {date})'))

        parts = [f'Rescan of `{image}` after the no-cache rebuild in {rebuild_link}.']
        if present is not None:
            if not survivors and not lingering:
                parts.append('The rebuild cleared every advisory it was dispatched for. '
                             + (f'Cleared by the rebuild: {", ".join(cleared)}.' if cleared else ''))
            else:
                parts.append(f'Cleared by the rebuild: {", ".join(cleared)}.' if cleared
                             else 'The rebuild cleared none of the advisories it was dispatched for.')
            if lingering:
                parts.append(f'Still in the image, but no longer a fixable HIGH/CRITICAL finding: '
                             f'{", ".join(lingering)}.')
        if survivors:
            parts.append('Still present after the rebuild:\n\n' + '\n'.join(
                f'- {v}: ' + (f'#{tracking[v]["number"]}' if tracking[v] else 'no issue yet') for v in survivors))
        parts.extend(paragraph for _, paragraph, _ in problems)
        if not problems and survivors:
            parts.append('The issues above track the survivors from here, so this cycle is done.')
        parts.extend(SCAN_RED_MARKER.format(marker) for marker, _, _ in problems)
        body = '\n\n'.join(part.strip() for part in parts)

        if not problems:
            return decide('close' if open_issue else 'none', body)
        if open_issue:
            news = any(marker not in markers for marker, _, _ in problems)
            return decide('comment', body) if news else decide('none', body)
        return decide('open', body, problems[0][2])

    if route == 'rebuild':
        retry = sorted(set(retry_ids) - set(ids))
        marker = 'ids=' + ','.join(sorted(set(ids) | set(retry)))
        parts = []
        if ids:
            parts.append(f'{run_link} found {_advisories(len(ids))} in `{image}`:\n\n{findings_table(report, ids)}')
        if retry:
            lead = 'Ubuntu has also published fixes, not yet in the image,' if ids else (
                f'{run_link} found that Ubuntu has published fixes, not yet in `{image}`,')
            parts.append(f'{lead} for advisories that open `cve` issues track:\n\n{findings_table(report, retry)}')
        parts += [DISPATCH_TEXT.get(dispatch, DISPATCH_FAILED_TEXT), LIFECYCLE_TEXT, SCAN_RED_MARKER.format(marker)]
        body = '\n\n'.join(parts)
        if open_issue:
            latest = next((m for m in reversed(markers) if m.startswith('ids=')), None)
            return decide('none', body) if latest == marker else decide('comment', body)
        title = (f'New fixable advisories in the main images (since {date})' if ids
                 else f'Ubuntu published fixes for open CVE issues; rebuilding the main images (since {date})')
        return decide('open', body, title)

    if route == 'none' and open_issue:
        if retry_ids:
            return decide('none', f'Nothing new as of {run_link}, but Ubuntu has published fixes for '
                                  f'{", ".join(sorted(retry_ids))} that this run could not rebuild for, so this '
                                  f'cycle stays open.')
        # A problem the current cycle reported (after its latest `ids=` marker) is still there.
        latest = max((n for n, m in enumerate(markers) if m.startswith('ids=')), default=-1)
        held = [m for m in markers[latest + 1:] if m.startswith(('rebuild-failed:', 'deploy-failed:'))]
        if held:
            return decide('none', f'Nothing new as of {run_link}, but this cycle still needs a maintainer '
                                  f'({", ".join(held)}), so it stays open.')
        body = (f'No new fixable HIGH/CRITICAL advisories remain in `{image}` as of {run_link}, '
                f'so this rebuild cycle is over.')
        return decide('close', body)
    return decide('none')


# Comment and issue authors as gh prints them: `gh issue list --json author` prefixes apps
# with "app/", comments do not.
BOT_LOGINS = frozenset({'github-actions', 'app/github-actions', 'github-actions[bot]'})
RETRY_MARKER = '<!-- cve-retry:{} -->'
RETRY_MARKER_RE = re.compile(r'<!-- cve-retry:(\S+?) -->')
# On the comments that close an issue as gone and that reopen it.
BACKLOG_MARKER = '<!-- cve-backlog:{} -->'
BACKLOG_MARKER_RE = re.compile(r'<!-- cve-backlog:(\S+?) -->')
# Debian package names and versions (Debian Policy 5.6.1 and 5.6.12). Anything else is
# dropped before it reaches the apt probe's shell.
DEB_NAME = re.compile(r'[a-z0-9][a-z0-9+.-]+')
DEB_VERSION = re.compile(r'[0-9][A-Za-z0-9.+~:-]*')
IMAGE_DIGEST = re.compile(r'sha256:[0-9a-f]{64}')
# update-cve-issues comments, then closes, seconds apart. A later close is someone else's.
GONE_CLOSE_WINDOW = datetime.timedelta(minutes=10)


def _login(author):
    return author.get('login') if isinstance(author, dict) else author


def _filed_by_this_workflow(issue):
    labels = {label.get('name') if isinstance(label, dict) else label for label in issue.get('labels') or []}
    return _login(issue.get('author')) in BOT_LOGINS and 'test' not in labels


def _bot_markers(issue, marker_re):
    """Markers in the comments this workflow wrote, oldest first: anyone can comment on a public repo."""
    return [m for m, _ in _bot_marker_comments(issue, marker_re)]


def _bot_marker_comments(issue, marker_re):
    """(marker, the comment's createdAt) for each marker in a comment this workflow wrote, oldest first."""
    return [(m, c.get('createdAt')) for c in issue.get('comments') or [] if _login(c.get('author')) in BOT_LOGINS
            for m in marker_re.findall(c.get('body') or '')]


def _time(value):
    """A GitHub API timestamp ("2026-11-06T09:10:02Z") as a datetime, or None."""
    try:
        return datetime.datetime.fromisoformat(value.replace('Z', '+00:00'))
    except (AttributeError, ValueError):
        return None


def _scanned_digest(report):
    """The digest of the image a Trivy report scanned (Metadata.RepoDigests), or None."""
    for ref in ((report or {}).get('Metadata') or {}).get('RepoDigests') or []:
        digest = ref.rpartition('@')[2]
        if IMAGE_DIGEST.fullmatch(digest):
            return digest
    return None


def _advisory_names(report):
    """Every ID and CVE/GHSA alias of every finding in a report, upper-cased."""
    return {a.upper() for result in (report or {}).get('Results') or []
            for vuln in result.get('Vulnerabilities') or []
            for a in [vuln.get('VulnerabilityID') or '', *(vuln.get('VendorIDs') or [])] if ADVISORY_ID.fullmatch(a)}


def backlog(issues, report, published=None, closed_issues=()):
    """Sort the `cve` issues this workflow filed by what report says about them. See the module doc.

    published is {(Debian package, fixed version): candidate} for each fixed version that
    today's Ubuntu archive has a candidate at least as new as, or None if the apt probe did
    not run.
    """
    results = (report or {}).get('Results') or []
    complete = any(r.get('Class') == 'os-pkgs' for r in results)
    findings = []  # (every ID and alias, upper-cased; is an Ubuntu package; the finding)
    for result in results:
        is_os = result.get('Class') == 'os-pkgs'
        for vuln in result.get('Vulnerabilities') or []:
            names = {a.upper() for a in [vuln.get('VulnerabilityID') or '', *(vuln.get('VendorIDs') or [])]
                     if ADVISORY_ID.fullmatch(a)}
            findings.append((names, is_os, vuln))

    gone, waiting, retry = [], [], []
    for issue in sorted(issues, key=lambda i: i['number']):
        ids = list(dict.fromkeys(ADVISORY_ID.findall(issue.get('title') or '')))
        if not ids or not _filed_by_this_workflow(issue):
            continue
        wanted = {i.upper() for i in ids}
        matched = [(is_os, vuln) for names, is_os, vuln in findings if names & wanted]
        if not matched:
            if complete:
                gone.append({'number': issue['number'], 'ids': ids})
            continue
        if not all(is_os for is_os, _ in matched):
            continue
        packages = {}
        for _, vuln in matched:
            name = vuln.get('PkgName') or ''
            fixed = [v.strip() for v in (vuln.get('FixedVersion') or '').split(',')]
            fixed = [v for v in fixed if DEB_VERSION.fullmatch(v)]
            if DEB_NAME.fullmatch(name) and fixed:
                versions = packages.setdefault(name, [])
                versions.extend(v for v in fixed if v not in versions)
        vuln_ids = sorted({vuln['VulnerabilityID'] for _, vuln in matched})
        waiting.append({'number': issue['number'], 'ids': ids, 'vuln_ids': vuln_ids,
                        'packages': [{'name': n, 'fixed': f} for n, f in sorted(packages.items())]})
        if published:
            tried = set(_bot_markers(issue, RETRY_MARKER_RE))
            markers = []
            for name, fixed in sorted(packages.items()):
                candidate = next((published[(name, f)] for f in fixed if (name, f) in published), None)
                if candidate:
                    markers.append(f'{name}={candidate}')
            untried = [m for m in markers if m not in tried]
            if untried:
                retry.append({'number': issue['number'], 'ids': ids, 'vuln_ids': vuln_ids, 'markers': untried})

    # A closed issue comes back only if its latest close was this workflow's close as gone,
    # which follows the gone comment within seconds: a fix PR's or a maintainer's close
    # stands, including one after reopening by hand, which leaves no marker.
    reopen = []
    for issue in sorted(closed_issues, key=lambda i: i['number']):
        ids = list(dict.fromkeys(ADVISORY_ID.findall(issue.get('title') or '')))
        if not ids or not _filed_by_this_workflow(issue):
            continue
        states = [(m, at) for m, at in _bot_marker_comments(issue, BACKLOG_MARKER_RE) if m in ('gone', 'reopened')]
        if not states or states[-1][0] != 'gone':
            continue
        gone_at, closed_at = _time(states[-1][1]), _time(issue.get('closedAt'))
        if gone_at is None or closed_at is None or not gone_at <= closed_at <= gone_at + GONE_CLOSE_WINDOW:
            continue
        wanted = {i.upper() for i in ids}
        if any(names & wanted for names, _, _ in findings):
            reopen.append({'number': issue['number'], 'ids': ids})
    return {'complete': complete, 'scanned_digest': _scanned_digest(report), 'gone': gone, 'waiting': waiting,
            'retry': retry, 'reopen': reopen, 'retry_ids': sorted({v for r in retry for v in r['vuln_ids']}),
            'retry_versions': sorted({m for r in retry for m in r['markers']})}


def backlog_actions(backlog, *, rebuild_state, post_rebuild_run_id='', pre_rebuild_versions=(), main_digest='',
                    repo_url='', run_id='', image=''):
    """The issue updates a backlog calls for: [{action: close|reopen|comment, number, body}]. See the module doc."""
    run_link = f'[run {run_id}]({repo_url}/actions/runs/{run_id})'
    actions = []
    # Gone or back is what the scanned image says. Once main-mega-amd64 points elsewhere (a
    # push since the scan, or since the rebuild a rescan pinned), the next daily scan decides.
    current = main_digest == backlog.get('scanned_digest')  # an unknown digest (None) never matches
    for entry in (backlog.get('gone') or []) if current else []:
        body = (f'{run_link} scanned `{image}` and no longer finds {", ".join(entry["ids"])}: whatever cleared '
                f'it (a newer upstream image, a no-cache rebuild, a fix PR or a `.trivyignore.yaml` entry) is in '
                f'main\'s images now, so this issue is done. If the advisory comes back, the scan reopens this '
                f'issue.\n\n' + BACKLOG_MARKER.format('gone'))
        actions.append({'action': 'close', 'number': entry['number'], 'body': body})
    for entry in (backlog.get('reopen') or []) if current else []:
        body = (f'{run_link} finds {", ".join(entry["ids"])} in `{image}` again. The scan closed this issue when '
                f'the advisory was gone, and it never files an advisory twice, so it reopens this one '
                f'instead.\n\n' + BACKLOG_MARKER.format('reopened'))
        actions.append({'action': 'reopen', 'number': entry['number'], 'body': body})
    if rebuild_state != 'rebuilt':
        # Before a rebuild the route dispatches one; after a failed one, apt never ran.
        return actions
    rebuild_link = f'[run {post_rebuild_run_id}]({repo_url}/actions/runs/{post_rebuild_run_id})'
    dispatched = set(pre_rebuild_versions)
    for entry in backlog.get('retry') or []:
        # Only the versions the dispatching scan saw published, so the rebuild could install
        # them: one Ubuntu published while it ran, or a rebuild run for other fixes, proves nothing.
        tried = [m for m in entry['markers'] if m in dispatched]
        if not tried:
            continue
        versions = ', '.join('`{}` {}'.format(*m.split('=', 1)) for m in tried)
        body = (f'The no-cache rebuild in {rebuild_link} was dispatched once Ubuntu\'s archive had {versions}, at '
                f'or above the fixed version, and ran `apt-get upgrade -y`, yet {run_link} still finds '
                f'{", ".join(entry["vuln_ids"])} in `{image}`. apt did not install the fix: it may be holding the '
                f'package back (an upgrade that needs a new dependency).\n\n**This needs a maintainer.** The '
                f'scan dispatches no further rebuild for these versions, only for a newer one Ubuntu publishes. '
                f'{RERUN_TEXT}\n\n' + '\n'.join(RETRY_MARKER.format(m) for m in tried))
        actions.append({'action': 'comment', 'number': entry['number'], 'body': body})
    return actions


def parse_published(text, out):
    """{(package, fixed version): candidate} from the apt probe's "<package> <fixed> <candidate>" lines."""
    published = {}
    for line in text.splitlines():
        fields = line.split()
        if (len(fields) == 3 and DEB_NAME.fullmatch(fields[0])
                and DEB_VERSION.fullmatch(fields[1]) and DEB_VERSION.fullmatch(fields[2])):
            published[(fields[0], fields[1])] = fields[2]
        elif fields:
            annotate(out, 'warning', f'Ignoring a malformed apt probe line: {line!r}')
    return published


def parse_bool(value):
    if value in ('true', 'false', ''):
        return value == 'true'
    raise argparse.ArgumentTypeError(f"expected 'true' or 'false', got {value!r}")


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
    p = commands.add_parser('route', help='where new IDs go: rebuild, triage or none')
    p.add_argument('--cve-ids', required=True, help='space-separated new IDs from new-ids')
    for flag in ('--test-mode', '--dry-run', '--skip-rebuild'):
        p.add_argument(flag, type=parse_bool, default=False, metavar='true|false')
    p.add_argument('--post-rebuild-run-id', default='', help='set when this run is a rebuild\'s rescan')
    p.add_argument('--retry-ids', default='', help='space-separated retry IDs from backlog')
    p.add_argument('--out', required=True, help='written with the route')
    p = commands.add_parser('scan-red', help='how to update the scan-red issue')
    p.add_argument('--route', required=True, choices=('rebuild', 'triage', 'none'))
    p.add_argument('--scan-result', default='success',
                   help='the scan job\'s result; only a rescan\'s failure is reported')
    p.add_argument('--cve-ids', required=True, help='space-separated new IDs from this run\'s scan')
    p.add_argument('--retry-ids', default='', help='space-separated retry IDs from backlog')
    p.add_argument('--pre-rebuild-ids', default='', help='the IDs the rebuild was dispatched for')
    p.add_argument('--post-rebuild-run-id', default='')
    p.add_argument('--rebuild-result', default='', help='result of the rebuild\'s build-mega-amd64 job')
    p.add_argument('--deploy-result', default='', help='result of the rebuild\'s deploy-to-quay job')
    p.add_argument('--dispatch', default='', help='request-rebuild output: dispatched or in-flight')
    p.add_argument('--file-issues-result', default='')
    p.add_argument('--open-issue-json', required=True,
                   help='the lowest-numbered open scan-red issue {number, body, comments: [{body}]}, or null')
    p.add_argument('--issues-json', required=True, help='every issue, open and closed: [{number, title}]')
    p.add_argument('--trivy-json', required=True, help='may be absent; the table then lists bare IDs')
    p.add_argument('--repo-url', required=True)
    p.add_argument('--run-id', required=True)
    p.add_argument('--image', required=True, help='the image ref this run scanned')
    p.add_argument('--date', required=True)
    p.add_argument('--out', required=True, help='written with the decision as JSON')
    p = commands.add_parser('backlog', help='what the scan says about the open cve issues')
    p.add_argument('--trivy-json', required=True)
    p.add_argument('--issues-json', required=True,
                   help='open cve issues: [{number, title, author: {login}, labels: [{name}], '
                        'comments: [{author: {login}, body}]}]')
    p.add_argument('--closed-issues-json', help='closed cve issues, in the same shape, to check for reopening')
    p.add_argument('--published',
                   help='the apt probe\'s "<package> <fixed> <candidate>" lines; omit if it did not run')
    p.add_argument('--waiting-pkgs-out', help='written with the "<package> <fixed version>" lines to probe')
    p.add_argument('--out', required=True, help='written with the backlog as JSON')
    p = commands.add_parser('backlog-actions', help='the issue updates a backlog calls for')
    p.add_argument('--backlog', required=True, help='backlog\'s JSON')
    p.add_argument('--rebuild-state', required=True, choices=('not-rebuilt', 'rebuilt', 'rebuild-failed'))
    p.add_argument('--post-rebuild-run-id', default='')
    p.add_argument('--pre-rebuild-versions', default='',
                   help='the "<package>=<version>" fixes the rebuild was dispatched for (backlog\'s retry_versions)')
    p.add_argument('--main-digest', default='',
                   help='what main-mega-amd64 points at now; empty if it could not be resolved')
    p.add_argument('--repo-url', required=True)
    p.add_argument('--run-id', required=True)
    p.add_argument('--image', required=True, help='the image ref this run scanned')
    p.add_argument('--out', required=True, help='written with the actions as JSON')
    args = parser.parse_args(argv)

    if args.command == 'names':
        return 0 if mentions(args.text, args.advisory_id) else 1
    if args.command == 'route':
        result = route(args.cve_ids.split(), args.test_mode, args.dry_run, args.skip_rebuild,
                       bool(args.post_rebuild_run_id), args.retry_ids.split())
        annotate(sys.stdout, 'notice', f'route={result} for new IDs: {args.cve_ids or "(none)"}'
                 + (f'; Ubuntu fixes published for: {args.retry_ids}' if args.retry_ids else ''))
        with open(args.out, 'w') as f:
            f.write(result + '\n')
        return 0
    if args.command == 'scan-red':
        with open(args.open_issue_json) as f:
            open_issue = json.load(f)
        with open(args.issues_json) as f:
            issues = json.load(f)
        try:
            with open(args.trivy_json) as f:
                trivy = json.load(f)
        except FileNotFoundError:
            annotate(sys.stdout, 'warning', f'{args.trivy_json} is missing; the scan-red table lists bare IDs')
            trivy = None
        decision = scan_red(
            route=args.route, cve_ids=args.cve_ids.split(), retry_ids=args.retry_ids.split(),
            pre_rebuild_ids=args.pre_rebuild_ids.split(),
            post_rebuild_run_id=args.post_rebuild_run_id, rebuild_result=args.rebuild_result,
            deploy_result=args.deploy_result, dispatch=args.dispatch,
            file_issues_result=args.file_issues_result, open_issue=open_issue, issues=issues, report=trivy,
            repo_url=args.repo_url, run_id=args.run_id, image=args.image, date=args.date,
            scan_result=args.scan_result)
        print(f'scan-red: {decision["action"]}' + (f' #{decision["number"]}' if decision['number'] else ''))
        with open(args.out, 'w') as f:
            json.dump(decision, f)
        return 0
    if args.command == 'backlog':
        with open(args.trivy_json) as f:
            trivy = json.load(f)
        with open(args.issues_json) as f:
            issues = json.load(f)
        published = None
        if args.published:
            with open(args.published) as f:
                published = parse_published(f.read(), sys.stdout)
        closed = []
        if args.closed_issues_json:
            with open(args.closed_issues_json) as f:
                closed = json.load(f)
        result = backlog(issues, trivy, published, closed)
        if not result['complete']:
            annotate(sys.stdout, 'warning', 'The Trivy report has no Ubuntu (os-pkgs) result, so no issue counts as gone')
        annotate(sys.stdout, 'notice', 'Open cve issues: '
                 f'gone {[e["number"] for e in result["gone"]]}, '
                 f'waiting on Ubuntu {[e["number"] for e in result["waiting"]]}, '
                 f'Ubuntu fix published {[e["number"] for e in result["retry"]]}, '
                 f'back again {[e["number"] for e in result["reopen"]]}')
        if args.waiting_pkgs_out:
            lines = dict.fromkeys(f'{p["name"]} {v}' for e in result['waiting'] for p in e['packages']
                                  for v in p['fixed'])
            with open(args.waiting_pkgs_out, 'w') as f:
                f.write(''.join(line + '\n' for line in lines))
        with open(args.out, 'w') as f:
            json.dump(result, f)
        return 0
    if args.command == 'backlog-actions':
        with open(args.backlog) as f:
            data = json.load(f)
        # A dispatched version counts only where it equals one of backlog's validated markers.
        result = backlog_actions(
            data, rebuild_state=args.rebuild_state, post_rebuild_run_id=args.post_rebuild_run_id,
            pre_rebuild_versions=args.pre_rebuild_versions.split(), main_digest=args.main_digest,
            repo_url=args.repo_url, run_id=args.run_id, image=args.image)
        held = [e['number'] for key in ('gone', 'reopen') for e in data.get(key) or []]
        if held and args.main_digest != data.get('scanned_digest'):
            annotate(sys.stdout, 'warning',
                     f'main-mega-amd64 is {args.main_digest or "unresolved"}, not the scanned '
                     f'{data.get("scanned_digest") or "(unknown digest)"}, so this run leaves '
                     f'{", ".join(f"#{n}" for n in held)} for the next daily scan to close or reopen')
        print('backlog actions: ' + (', '.join(f'{a["action"]} #{a["number"]}' for a in result) or 'none'))
        with open(args.out, 'w') as f:
            json.dump(result, f)
        return 0
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
