"""Unit tests for cve_pipeline.py (no network access).

Run with: python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
"""

import io
import itertools
import json
import os
import tempfile
import unittest
from unittest import mock

import cve_pipeline as cp


def finding(vuln_id, severity='HIGH', fixed='1.2.3', vendor_ids=None):
    vuln = {'VulnerabilityID': vuln_id, 'PkgName': 'pkg', 'Severity': severity, 'FixedVersion': fixed}
    if vendor_ids is not None:
        vuln['VendorIDs'] = vendor_ids
    return vuln


def report(*vulns):
    return {'Results': [{'Target': 'image', 'Vulnerabilities': list(vulns)}]}


def issue(number, title):
    return {'number': number, 'title': title}


# Modeled on the five issues PR #1137 closed.
ISSUES_1137 = [
    issue(1132, '[CVE-2026-112403] pillow: heap buffer overflow in JPEG draft-scale decoding (potential RCE)'),
    issue(1133, '[CVE-2026-112404] pillow: out-of-bounds read in TIFF tile/RGBA decoding'),
    issue(1134, '[CVE-2026-112896] pillow: out-of-bounds read in TIFF encoder when re-saving a crafted ColorMap'),
    issue(1135, '[CVE-2026-112897] pillow: uninitialized heap memory disclosure in short TIFF strip/tile decodes'),
    issue(1136, '[GHSA-6q48-6c7w-9gw2] certifi: bundled emSign root CA chains to compromised intermediate CA keys'),
]
BODY_1137 = '''## Summary

Bumps Pillow to 13.0.0 and floors certifi.

## Closes

Closes #1132
Closes #1133
Closes #1134
Closes #1135
Closes #1136

## Per-CVE rationale
'''
GOOD_TITLE = ('Bump pillow to 13.0.0, floor certifi (CVE-2026-112403, CVE-2026-112404, '
              'CVE-2026-112896, CVE-2026-112897, GHSA-6q48-6c7w-9gw2)')
FALLBACK_1137 = ('Auto-fix: pillow CVE-2026-112403, CVE-2026-112404, CVE-2026-112896, CVE-2026-112897; '
                 'certifi GHSA-6q48-6c7w-9gw2')


class MentionsTest(unittest.TestCase):

    def test_bracketed_id_is_mentioned(self):
        self.assertTrue(cp.mentions('[CVE-2026-1001] foo: bar', 'CVE-2026-1001'))

    def test_longer_id_does_not_mention_its_prefix(self):
        self.assertFalse(cp.mentions('[CVE-2026-10011] foo: bar', 'CVE-2026-1001'))

    def test_id_preceded_by_a_letter_is_not_mentioned(self):
        self.assertFalse(cp.mentions('[NOTCVE-2026-1001] foo: bar', 'CVE-2026-1001'))

    def test_match_ignores_case(self):
        self.assertTrue(cp.mentions('[ghsa-6q48-6c7w-9gw2] certifi: x', 'GHSA-6q48-6c7w-9gw2'))


class NewAdvisoriesTest(unittest.TestCase):

    def new(self, scan, issues):
        out = io.StringIO()
        return cp.new_advisories(scan, issues, out), out.getvalue()

    def test_untracked_finding_is_new(self):
        new, _ = self.new(report(finding('CVE-2026-1001')), [])
        self.assertEqual(new, ['CVE-2026-1001'])

    def test_finding_named_in_an_issue_title_is_not_new(self):
        new, out = self.new(report(finding('CVE-2026-1001')), [issue(10, '[CVE-2026-1001] foo: bar')])
        self.assertEqual(new, [])
        self.assertIn('#10', out)

    def test_only_fixable_high_and_critical_findings_count(self):
        scan = report(
            finding('CVE-2026-1001', severity='MEDIUM'),
            finding('CVE-2026-1002', fixed=''),
            {'VulnerabilityID': 'CVE-2026-1003', 'Severity': 'HIGH'},
            finding('CVE-2026-1004', severity='CRITICAL'),
        )
        new, _ = self.new(scan, [])
        self.assertEqual(new, ['CVE-2026-1004'])

    def test_id_that_prefixes_a_tracked_id_is_still_new(self):
        new, _ = self.new(report(finding('CVE-2026-1123')), [issue(10, '[CVE-2026-11234] foo: bar')])
        self.assertEqual(new, ['CVE-2026-1123'])

    def test_advisory_rekeyed_from_ghsa_to_cve_is_not_filed_again(self):
        # Trivy reports a GitHub advisory under its CVE once one is assigned, and moves
        # the GHSA ID to VendorIDs.
        scan = report(finding('CVE-2026-2002', vendor_ids=['GHSA-6q48-6c7w-9gw2']))
        new, out = self.new(scan, [ISSUES_1137[4]])
        self.assertEqual(new, [])
        for expected in ('CVE-2026-2002', 'GHSA-6q48-6c7w-9gw2', '#1136', '.trivyignore.yaml'):
            self.assertIn(expected, out)

    def test_aliases_are_pooled_across_findings_of_the_same_id(self):
        scan = {'Results': [
            {'Target': 'a', 'Vulnerabilities': [finding('CVE-2026-2002')]},
            {'Target': 'b', 'Vulnerabilities': [finding('CVE-2026-2002', vendor_ids=['GHSA-6q48-6c7w-9gw2'])]},
        ]}
        new, _ = self.new(scan, [ISSUES_1137[4]])
        self.assertEqual(new, [])

    def test_only_cve_and_ghsa_vendor_ids_count_as_aliases(self):
        # Distro IDs such as Debian's DSA-* are not how issues are titled.
        scan = report(finding('CVE-2026-3003', vendor_ids=['DSA-5555-1']))
        new, _ = self.new(scan, [issue(30, '[DSA-5555-1] tracked under a distro ID')])
        self.assertEqual(new, ['CVE-2026-3003'])

    def test_new_ids_are_sorted_and_unique(self):
        scan = report(finding('CVE-2026-1002'), finding('CVE-2026-1001'), finding('CVE-2026-1002'))
        new, _ = self.new(scan, [])
        self.assertEqual(new, ['CVE-2026-1001', 'CVE-2026-1002'])

    def test_results_without_vulnerabilities(self):
        scan = {'Results': [{'Target': 'a', 'Vulnerabilities': None}, {'Target': 'b'}]}
        new, _ = self.new(scan, [])
        self.assertEqual(new, [])

    def test_warning_is_a_workflow_annotation_under_actions(self):
        scan = report(finding('CVE-2026-2002', vendor_ids=['GHSA-6q48-6c7w-9gw2']))
        with mock.patch.dict(os.environ, {'GITHUB_ACTIONS': 'true'}):
            _, out = self.new(scan, [ISSUES_1137[4]])
        self.assertIn('::warning::', out)


class PrTitleTest(unittest.TestCase):

    def title(self, agent_title, body=BODY_1137, issues=ISSUES_1137):
        out = io.StringIO()
        return cp.pr_title(agent_title, body, issues, out), out.getvalue()

    def test_agent_title_naming_every_closed_advisory_is_used(self):
        title, out = self.title(GOOD_TITLE + '\n')
        self.assertEqual(title, GOOD_TITLE)
        self.assertNotIn('WARNING', out.upper())

    def test_no_title_file_falls_back(self):
        title, out = self.title(None)
        self.assertEqual(title, FALLBACK_1137)
        self.assertIn('TITLE.txt', out)

    def test_title_missing_an_advisory_falls_back(self):
        title, out = self.title('Bump pillow to 13.0.0 (CVE-2026-112403, CVE-2026-112404, '
                                'CVE-2026-112896, CVE-2026-112897)')
        self.assertEqual(title, FALLBACK_1137)
        self.assertIn('GHSA-6q48-6c7w-9gw2', out)

    def test_longer_id_in_title_does_not_name_its_prefix(self):
        title, _ = self.title('Bump foo (CVE-2026-11234)', body='Closes #5',
                              issues=[issue(5, '[CVE-2026-1123] foo: bar')])
        self.assertEqual(title, 'Auto-fix: foo CVE-2026-1123')

    def test_multiline_title_falls_back(self):
        title, _ = self.title(GOOD_TITLE + '\n\nAlso fixes the build.')
        self.assertEqual(title, FALLBACK_1137)

    def test_blank_title_falls_back(self):
        title, _ = self.title('  \n')
        self.assertEqual(title, FALLBACK_1137)

    def test_overlong_title_falls_back(self):
        title, out = self.title(GOOD_TITLE + ' ' + 'x' * cp.MAX_TITLE_LENGTH)
        self.assertEqual(title, FALLBACK_1137)
        self.assertIn(str(cp.MAX_TITLE_LENGTH), out)

    def test_mentions_links_and_markup_fall_back(self):
        for extra in (' @someone', ' `pillow`', ' #1137', ' see https://example.com', ' [x](y)', ' *now*'):
            with self.subTest(extra=extra):
                title, _ = self.title(GOOD_TITLE + extra)
                self.assertEqual(title, FALLBACK_1137)

    def test_version_specifiers_are_allowed(self):
        agent = ("Pin certifi>=2026.10.19 and pillow~=13.0, drop pip's vendored copy_ "
                 "(CVE-2026-112403, CVE-2026-112404, CVE-2026-112896, CVE-2026-112897, GHSA-6q48-6c7w-9gw2)")
        title, _ = self.title(agent)
        self.assertEqual(title, agent)

    def test_any_github_closing_keyword_counts(self):
        title, _ = self.title('Bump pillow (CVE-2026-112403)', body='Fixes #1132\nresolves #1136\n')
        self.assertEqual(title, 'Auto-fix: pillow CVE-2026-112403; certifi GHSA-6q48-6c7w-9gw2')

    def test_closed_issue_that_is_not_an_open_cve_issue_is_reported_and_skipped(self):
        title, out = self.title('Bump pillow (CVE-2026-112403)', body='Closes #1132\nCloses #999\n')
        self.assertEqual(title, 'Bump pillow (CVE-2026-112403)')
        self.assertIn('#999', out)

    def test_body_closing_nothing_gets_the_generic_title(self):
        for agent in (None, GOOD_TITLE):
            with self.subTest(agent=agent):
                title, out = self.title(agent, body='## Summary\n\nNo closing lines.\n')
                self.assertEqual(title, cp.GENERIC_TITLE)
                self.assertIn('closes no', out)

    def test_fallback_groups_by_package_in_closing_order(self):
        title, _ = self.title(None, body='Closes #1136\nCloses #1132\nCloses #1133\n')
        self.assertEqual(title, 'Auto-fix: certifi GHSA-6q48-6c7w-9gw2; pillow CVE-2026-112403, CVE-2026-112404')

    def test_fallback_shortens_when_the_ids_do_not_fit(self):
        issues = [issue(100 + n, f'[CVE-2026-1000{n:02d}] pillow: bug {n}') for n in range(15)]
        body = '\n'.join(f'Closes #{i["number"]}' for i in issues)
        title, _ = self.title(None, body=body, issues=issues)
        self.assertEqual(title, 'Auto-fix: pillow (15 advisories)')

    def test_fallback_keeps_the_ids_when_only_the_package_names_do_not_fit(self):
        # IDs outrank package names: release notes and `git log --grep CVE-` both need them.
        issues = [issue(200 + n, f'[CVE-2026-1000{n:02d}] org.example.very.long.group.name{n}:artifact-{n}: bug')
                  for n in range(7)]
        body = '\n'.join(f'Closes #{i["number"]}' for i in issues)
        title, _ = self.title(None, body=body, issues=issues)
        self.assertEqual(title, 'Auto-fix: ' + ', '.join(f'CVE-2026-1000{n:02d}' for n in range(7)))

    def test_fallback_keeps_maven_coordinates_whole(self):
        issues = [issue(7, '[CVE-2025-68161] org.apache.logging.log4j:log4j-core: socket appender skips TLS hostname check')]
        title, _ = self.title(None, body='Closes #7', issues=issues)
        self.assertEqual(title, 'Auto-fix: org.apache.logging.log4j:log4j-core CVE-2025-68161')

    def test_fallback_omits_a_package_it_cannot_parse_or_trust(self):
        for issue_title in ('CVE-2026-1001 in some package', '[CVE-2026-1001] `pkg`@evil: bug'):
            with self.subTest(issue_title=issue_title):
                title, _ = self.title(None, body='Closes #8', issues=[issue(8, issue_title)])
                self.assertEqual(title, 'Auto-fix: CVE-2026-1001')


class RouteTest(unittest.TestCase):
    """Where the scan's new IDs go: a no-cache rebuild first, Claude triage, or nowhere."""

    FLAGS = ('test_mode', 'dry_run', 'skip_rebuild', 'post_rebuild')

    def test_new_ids_on_a_plain_run_go_to_a_rebuild_first(self):
        self.assertEqual(cp.route(['CVE-2026-1001']), 'rebuild')

    def test_no_new_ids_go_nowhere(self):
        for flag in ('dry_run', 'skip_rebuild', 'post_rebuild'):
            with self.subTest(flag=flag):
                self.assertEqual(cp.route([], **{flag: True}), 'none')
        self.assertEqual(cp.route([]), 'none')

    def test_test_mode_triages_directly(self):
        self.assertEqual(cp.route(['CVE-2026-1001'], test_mode=True), 'triage')

    def test_dry_run_and_skip_rebuild_triage_directly(self):
        # A dry run must not set off a real rebuild that republishes main's images.
        for flag in ('dry_run', 'skip_rebuild'):
            with self.subTest(flag=flag):
                self.assertEqual(cp.route(['CVE-2026-1001'], **{flag: True}), 'triage')

    def test_survivors_of_a_rebuild_are_triaged(self):
        self.assertEqual(cp.route(['CVE-2026-1001'], post_rebuild=True), 'triage')

    def test_a_post_rebuild_run_never_routes_to_another_rebuild(self):
        # Loop safety: the rebuild's rescan must not dispatch a rebuild of its own.
        for values in itertools.product((False, True), repeat=len(self.FLAGS)):
            flags = dict(zip(self.FLAGS, values))
            if not flags['post_rebuild']:
                continue
            for ids in ([], ['CVE-2026-1001']):
                with self.subTest(ids=ids, **flags):
                    self.assertNotEqual(cp.route(ids, **flags), 'rebuild')


def os_finding(vuln_id, pkg, installed, fixed, severity='HIGH'):
    return {'VulnerabilityID': vuln_id, 'PkgName': pkg, 'InstalledVersion': installed,
            'FixedVersion': fixed, 'Severity': severity}


APT_REPORT = {'Results': [
    {'Target': 'ghcr.io/broadinstitute/viral-ngs:main-mega-amd64 (ubuntu 26.04)', 'Class': 'os-pkgs',
     'Type': 'ubuntu', 'Vulnerabilities': [
         os_finding('CVE-2026-109518', 'apt', '3.2.0', '3.2.1ubuntu0.1'),
         os_finding('CVE-2026-109518', 'libapt-pkg7.0', '3.2.0', '3.2.1ubuntu0.1'),
         os_finding('CVE-2026-5555', 'libpcre2-8-0', '10.46-1build1', '10.46-1ubuntu0.1', severity='MEDIUM'),
     ]},
    {'Target': 'opt/conda/conda-meta/pillow.json', 'Class': 'lang-pkgs', 'Type': 'conda-pkg',
     'Vulnerabilities': [os_finding('CVE-2026-7777', 'pillow', '12.3.0', '13.0.0', severity='CRITICAL')]},
]}
REPO = 'https://github.com/broadinstitute/viral-ngs'


def open_issue(number, body='', comments=()):
    return {'number': number, 'body': body, 'comments': [{'body': c} for c in comments]}


class ScanRedTest(unittest.TestCase):
    """The scan-red issue tracks one rebuild cycle: opened when the daily scan dispatches a
    no-cache rebuild, closed by that rebuild's rescan, left open when something needs a person."""

    def decide(self, **kwargs):
        args = dict(route='none', cve_ids=[], report=APT_REPORT, repo_url=REPO, run_id='900',
                    image='ghcr.io/broadinstitute/viral-ngs:main-mega-amd64', date='2026-11-05')
        args.update(kwargs)
        return cp.scan_red(**args)

    def rescan(self, **kwargs):
        args = dict(route='triage', post_rebuild_run_id='800', rebuild_result='success',
                    pre_rebuild_ids=['CVE-2026-109518', 'CVE-2026-7777'], file_issues_result='skipped')
        args.update(kwargs)
        return self.decide(**args)

    # --- the daily scan dispatches a rebuild ---

    def test_a_dispatched_rebuild_opens_an_issue(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-109518', 'CVE-2026-7777'], dispatch='dispatched')
        self.assertEqual(d['action'], 'open')
        self.assertIsNone(d['number'])
        self.assertIn('2026-11-05', d['title'])
        self.assertIn(f'{REPO}/actions/runs/900', d['body'])
        self.assertIn('force_rebuild', d['body'])
        self.assertIn('<!-- scan-red:ids=CVE-2026-109518,CVE-2026-7777 -->', d['body'])

    def test_the_issue_tables_each_new_finding(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-109518', 'CVE-2026-7777'], dispatch='dispatched')
        self.assertRegex(d['body'], r'\| CVE-2026-109518 \| HIGH \| apt, libapt-pkg7\.0 \| ubuntu \| 3\.2\.0 → 3\.2\.1ubuntu0\.1 \|')
        self.assertRegex(d['body'], r'\| CVE-2026-7777 \| CRITICAL \| pillow \| conda-pkg \| 12\.3\.0 → 13\.0\.0 \|')
        # Only this run's new IDs: not the MEDIUM finding, which is never triaged.
        self.assertNotIn('CVE-2026-5555', d['body'])

    def test_table_cells_cannot_break_the_table(self):
        scan = report(os_finding('CVE-2026-1001', 'evil|pkg', '1', '2'))
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001'], dispatch='dispatched', report=scan)
        self.assertIn('evil\\|pkg', d['body'])

    def test_a_missing_trivy_report_still_lists_the_ids(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001'], dispatch='dispatched', report=None)
        self.assertEqual(d['action'], 'open')
        self.assertIn('CVE-2026-1001', d['body'])

    def test_a_rebuild_already_in_flight_is_said_so(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001'], dispatch='in-flight')
        self.assertEqual(d['action'], 'open')
        self.assertIn('already in flight', d['body'])

    def test_a_failed_dispatch_is_said_so(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001'], dispatch='')
        self.assertEqual(d['action'], 'open')
        self.assertIn('could not be dispatched', d['body'])

    def test_a_stuck_cycle_gets_a_comment_when_the_ids_change(self):
        stuck = open_issue(1150, body='<!-- scan-red:ids=CVE-2026-1001 -->')
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001', 'CVE-2026-1002'], dispatch='dispatched',
                        open_issue=stuck)
        self.assertEqual((d['action'], d['number']), ('comment', 1150))
        self.assertIn('<!-- scan-red:ids=CVE-2026-1001,CVE-2026-1002 -->', d['body'])

    def test_a_stuck_cycle_is_not_commented_on_daily_for_the_same_ids(self):
        stuck = open_issue(1150, body='<!-- scan-red:ids=CVE-2026-1001 -->',
                           comments=['<!-- scan-red:ids=CVE-2026-1001,CVE-2026-1002 -->'])
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1002', 'CVE-2026-1001'], dispatch='dispatched',
                        open_issue=stuck)
        self.assertEqual(d['action'], 'none')

    def test_dedup_compares_against_the_latest_marker(self):
        # The ID set went back to an earlier one: that is a change from the latest note.
        stuck = open_issue(1150, body='<!-- scan-red:ids=CVE-2026-1001 -->',
                           comments=['<!-- scan-red:ids=CVE-2026-1001,CVE-2026-1002 -->'])
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-1001'], dispatch='dispatched', open_issue=stuck)
        self.assertEqual(d['action'], 'comment')

    # --- the rebuild's rescan ---

    def test_a_rebuild_that_cleared_everything_closes_the_issue(self):
        d = self.rescan(cve_ids=[], open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        for expected in ('CVE-2026-109518', 'CVE-2026-7777', f'{REPO}/actions/runs/800'):
            self.assertIn(expected, d['body'])

    def test_a_rebuild_nobody_was_told_about_closes_nothing(self):
        # A maintainer's own force_rebuild dispatch: no issue is open, and none is needed.
        d = self.rescan(cve_ids=[])
        self.assertEqual(d['action'], 'none')
        self.assertIn('CVE-2026-109518', d['summary'])

    def test_survivors_filed_as_cve_issues_close_the_cycle_with_links(self):
        filed = [issue(1151, '[CVE-2026-7777] pillow: out-of-bounds read')]
        d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result='success', cve_issues=filed,
                        open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('#1151', d['body'])
        self.assertRegex(d['body'], r'(?s)[Cc]leared.*CVE-2026-109518')

    def test_survivors_that_were_not_filed_keep_the_issue_open(self):
        for result in ('failure', 'skipped', 'cancelled'):
            with self.subTest(file_issues_result=result):
                d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result=result, open_issue=open_issue(1150))
                self.assertEqual((d['action'], d['number']), ('comment', 1150))
                self.assertIn('CVE-2026-7777', d['body'])

    def test_survivors_that_were_not_filed_open_an_issue_if_none_is_open(self):
        # The rescan was started by github-actions[bot], so a failed run notifies nobody.
        d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result='failure')
        self.assertEqual(d['action'], 'open')

    def test_a_failed_rebuild_needs_a_maintainer(self):
        for result in ('failure', 'cancelled', 'skipped', ''):
            with self.subTest(rebuild_result=result):
                d = self.rescan(cve_ids=['CVE-2026-109518'], rebuild_result=result, open_issue=open_issue(1150))
                self.assertEqual((d['action'], d['number']), ('comment', 1150))
                self.assertIn('maintainer', d['body'])
                self.assertIn('<!-- scan-red:rebuild-failed:800 -->', d['body'])

    def test_a_failed_rebuild_is_reported_once(self):
        noted = open_issue(1150, comments=['<!-- scan-red:rebuild-failed:800 -->'])
        d = self.rescan(cve_ids=['CVE-2026-109518'], rebuild_result='failure', open_issue=noted)
        self.assertEqual(d['action'], 'none')

    def test_a_failed_rebuild_opens_an_issue_if_none_is_open(self):
        d = self.rescan(cve_ids=['CVE-2026-109518'], rebuild_result='failure')
        self.assertEqual(d['action'], 'open')
        self.assertIn('maintainer', d['body'])

    # --- the daily scan with nothing new ---

    def test_a_daily_scan_with_nothing_new_closes_a_leftover_issue(self):
        d = self.decide(route='none', open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('No new fixable', d['body'])

    def test_a_daily_scan_with_nothing_new_and_nothing_open_does_nothing(self):
        self.assertEqual(self.decide(route='none')['action'], 'none')

    def test_a_direct_triage_leaves_the_issue_alone(self):
        # skip_rebuild: a maintainer sent these straight to Claude; no rebuild cycle exists.
        d = self.decide(route='triage', cve_ids=['CVE-2026-1001'], open_issue=open_issue(1150))
        self.assertEqual(d['action'], 'none')


class CliTest(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def path(self, name, content=None):
        p = os.path.join(self.dir.name, name)
        if content is not None:
            with open(p, 'w') as f:
                f.write(content if isinstance(content, str) else json.dumps(content))
        return p

    def read(self, name):
        with open(self.path(name)) as f:
            return f.read()

    def run_main(self, argv):
        with mock.patch('sys.stdout', io.StringIO()):
            return cp.main(argv)

    def test_new_ids(self):
        rc = self.run_main([
            'new-ids',
            '--trivy-json', self.path('trivy.json', report(finding('CVE-2026-1002'), finding('CVE-2026-1001'))),
            '--issues-json', self.path('issues.json', []),
            '--out', self.path('new-ids.txt'),
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(self.read('new-ids.txt'), 'CVE-2026-1001 CVE-2026-1002\n')

    def test_pr_title_uses_the_agent_title(self):
        rc = self.run_main([
            'pr-title',
            '--title-file', self.path('TITLE.txt', GOOD_TITLE + '\n'),
            '--pr-body', self.path('PROCEED.md', BODY_1137),
            '--issues-json', self.path('issues.json', ISSUES_1137),
            '--out', self.path('pr-title.txt'),
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(self.read('pr-title.txt'), GOOD_TITLE + '\n')

    def test_names_matches_whole_ids_only(self):
        # The file-issues job's title guard: a title dedup won't recognize is refused.
        cases = [
            ('CVE-2026-1123', '[CVE-2026-1123] foo: bar', 0),
            ('CVE-2026-1123', '[CVE-2026-11234] foo: bar', 1),
            ('GHSA-6q48-6c7w-9gw2', '[ghsa-6q48-6c7w-9gw2] certifi: bar', 0),
        ]
        for advisory_id, text, rc in cases:
            with self.subTest(advisory_id=advisory_id, text=text):
                self.assertEqual(self.run_main(['names', '--id', advisory_id, '--text', text]), rc)

    def test_pr_title_without_a_title_file(self):
        rc = self.run_main([
            'pr-title',
            '--title-file', self.path('TITLE.txt'),
            '--pr-body', self.path('PROCEED.md', BODY_1137),
            '--issues-json', self.path('issues.json', ISSUES_1137),
            '--out', self.path('pr-title.txt'),
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(self.read('pr-title.txt'), FALLBACK_1137 + '\n')

    def test_route(self):
        # Workflow booleans arrive as the strings 'true' and 'false'.
        cases = [
            (['--cve-ids', 'CVE-2026-1001'], 'rebuild'),
            (['--cve-ids', ''], 'none'),
            (['--cve-ids', 'CVE-2026-1001', '--test-mode', 'true'], 'triage'),
            (['--cve-ids', 'CVE-2026-1001', '--dry-run', 'true', '--skip-rebuild', 'false'], 'triage'),
            (['--cve-ids', 'CVE-2026-1001 CVE-2026-1002', '--post-rebuild-run-id', '800'], 'triage'),
        ]
        for args, expected in cases:
            with self.subTest(args=args):
                self.assertEqual(self.run_main(['route', *args, '--out', self.path('route.txt')]), 0)
                self.assertEqual(self.read('route.txt'), expected + '\n')

    def test_scan_red(self):
        rc = self.run_main([
            'scan-red',
            '--route', 'triage',
            '--cve-ids', '',
            '--pre-rebuild-ids', 'CVE-2026-109518 CVE-2026-7777',
            '--post-rebuild-run-id', '800',
            '--rebuild-result', 'success',
            '--dispatch', '',
            '--file-issues-result', 'skipped',
            '--open-issue-json', self.path('open.json', open_issue(1150)),
            '--cve-issues-json', self.path('cve-issues.json', []),
            '--trivy-json', self.path('trivy.json', APT_REPORT),
            '--repo-url', REPO,
            '--run-id', '900',
            '--image', 'ghcr.io/broadinstitute/viral-ngs@sha256:abc',
            '--date', '2026-11-05',
            '--out', self.path('decision.json'),
        ])
        self.assertEqual(rc, 0)
        decision = json.loads(self.read('decision.json'))
        self.assertEqual((decision['action'], decision['number']), ('close', 1150))

    def test_scan_red_with_no_open_issue_and_no_trivy_report(self):
        rc = self.run_main([
            'scan-red',
            '--route', 'rebuild',
            '--cve-ids', 'CVE-2026-1001',
            '--dispatch', 'dispatched',
            '--open-issue-json', self.path('open.json', 'null'),
            '--cve-issues-json', self.path('cve-issues.json', []),
            '--trivy-json', self.path('missing-trivy.json'),
            '--repo-url', REPO,
            '--run-id', '900',
            '--image', 'ghcr.io/broadinstitute/viral-ngs:main-mega-amd64',
            '--date', '2026-11-05',
            '--out', self.path('decision.json'),
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(self.read('decision.json'))['action'], 'open')


if __name__ == '__main__':
    unittest.main()
