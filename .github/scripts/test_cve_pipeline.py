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
            for ids, retry_ids in itertools.product(([], ['CVE-2026-1001']), ([], ['CVE-2026-109518'])):
                with self.subTest(ids=ids, retry_ids=retry_ids, **flags):
                    self.assertNotEqual(cp.route(ids, retry_ids=retry_ids, **flags), 'rebuild')

    def test_a_published_ubuntu_fix_for_an_open_issue_goes_to_a_rebuild(self):
        # Nothing new, but Ubuntu now has the fix an open os-package issue is waiting for.
        self.assertEqual(cp.route([], retry_ids=['CVE-2026-109518']), 'rebuild')
        self.assertEqual(cp.route(['CVE-2026-1001'], retry_ids=['CVE-2026-109518']), 'rebuild')

    def test_a_retry_is_never_triaged(self):
        # Its advisory already has an issue: when no rebuild may run, there is nothing to do.
        for flag in ('dry_run', 'skip_rebuild', 'post_rebuild'):
            with self.subTest(flag=flag):
                self.assertEqual(cp.route([], retry_ids=['CVE-2026-109518'], **{flag: True}), 'none')


def os_finding(vuln_id, pkg, installed, fixed, severity='HIGH'):
    return {'VulnerabilityID': vuln_id, 'PkgName': pkg, 'InstalledVersion': installed,
            'FixedVersion': fixed, 'Severity': severity}


# What main-mega-amd64 pointed at when APT_REPORT was scanned.
MAIN_DIGEST = 'sha256:' + 'a' * 64
APT_REPORT = {'Metadata': {'RepoDigests': [f'ghcr.io/broadinstitute/viral-ngs@{MAIN_DIGEST}']}, 'Results': [
    {'Target': 'ghcr.io/broadinstitute/viral-ngs:main-mega-amd64 (ubuntu 26.04)', 'Class': 'os-pkgs',
     'Type': 'ubuntu', 'Vulnerabilities': [
         os_finding('CVE-2026-109518', 'apt', '3.2.0', '3.2.1ubuntu0.1'),
         os_finding('CVE-2026-109518', 'libapt-pkg7.0', '3.2.0', '3.2.1ubuntu0.1'),
         os_finding('CVE-2026-5555', 'libpcre2-8-0', '10.46-1build1', '10.46-1ubuntu0.1', severity='MEDIUM'),
     ]},
    {'Target': 'opt/conda/conda-meta/pillow.json', 'Class': 'lang-pkgs', 'Type': 'conda-pkg',
     'Vulnerabilities': [os_finding('CVE-2026-7777', 'pillow', '12.3.0', '13.0.0', severity='CRITICAL')]},
]}
# The rescan of a rebuild that cleared both: only the MEDIUM finding, never triaged, is left.
REBUILT_CLEAN = {'Results': [
    {'Target': 'ghcr.io/broadinstitute/viral-ngs@sha256:abc (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu',
     'Vulnerabilities': [os_finding('CVE-2026-5555', 'libpcre2-8-0', '10.46-1build1', '10.46-1ubuntu0.1',
                                    severity='MEDIUM')]},
]}
# The rescan of a rebuild that cleared the apt finding but not the conda one.
REBUILT_PILLOW = {'Results': [
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
        args = dict(route='triage', post_rebuild_run_id='800', rebuild_result='success', deploy_result='success',
                    pre_rebuild_ids=['CVE-2026-109518', 'CVE-2026-7777'], file_issues_result='skipped',
                    report=REBUILT_CLEAN)
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
        d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result='success', issues=filed,
                        report=REBUILT_PILLOW, open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('#1151', d['body'])
        self.assertIn('Cleared by the rebuild: CVE-2026-109518.', d['body'])

    def test_survivors_that_were_not_filed_keep_the_issue_open(self):
        for result in ('failure', 'skipped', 'cancelled'):
            with self.subTest(file_issues_result=result):
                d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result=result, report=REBUILT_PILLOW,
                                open_issue=open_issue(1150))
                self.assertEqual((d['action'], d['number']), ('comment', 1150))
                self.assertIn('CVE-2026-7777', d['body'])

    def test_survivors_that_were_not_filed_open_an_issue_if_none_is_open(self):
        # The rescan was started by github-actions[bot], so a failed run notifies nobody.
        d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result='failure', report=REBUILT_PILLOW)
        self.assertEqual(d['action'], 'open')

    def test_a_partly_filed_batch_keeps_the_issue_open(self):
        # file-issues succeeds when triage stopped mid-batch: it files only the reports it got.
        scan = report(finding('CVE-2026-7777'), finding('CVE-2026-8888'))
        filed = [issue(1151, '[CVE-2026-7777] pillow: out-of-bounds read')]
        d = self.rescan(cve_ids=['CVE-2026-7777', 'CVE-2026-8888'], file_issues_result='success', issues=filed,
                        report=scan, open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('comment', 1150))
        self.assertIn('#1151', d['body'])
        self.assertRegex(d['body'], r'CVE-2026-8888: no issue')

    def test_clearance_comes_from_the_rescan_not_from_the_new_ids(self):
        # CVE-2026-7777 got an issue while the rebuild ran, so it is no longer "new", but the
        # rebuilt image still has it: it survived, it was not cleared.
        filed = [issue(1151, '[CVE-2026-7777] pillow: out-of-bounds read')]
        d = self.rescan(cve_ids=[], issues=filed, report=REBUILT_PILLOW, open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('Cleared by the rebuild: CVE-2026-109518.', d['body'])
        self.assertIn('CVE-2026-7777: #1151', d['body'])
        self.assertNotIn('cleared every advisory', d['body'])

    def test_a_survivor_reported_under_a_new_alias_is_matched_to_its_issue(self):
        # Dispatched for the GHSA; the rebuilt image's report lists it as the CVE it was
        # since assigned, with the GHSA in VendorIDs. Its issue is titled with the GHSA.
        scan = report(finding('CVE-2026-2002', vendor_ids=['GHSA-6q48-6c7w-9gw2']))
        d = self.rescan(pre_rebuild_ids=['GHSA-6q48-6c7w-9gw2'], cve_ids=[], issues=[ISSUES_1137[4]],
                        report=scan, open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('CVE-2026-2002: #1136', d['body'])
        self.assertNotIn('Cleared by the rebuild: GHSA', d['body'])

    def test_a_missing_rescan_report_claims_nothing_cleared(self):
        d = self.rescan(cve_ids=[], report=None, open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('comment', 1150))
        self.assertNotIn('cleared every advisory', d['body'])
        self.assertIn('Trivy report', d['body'])
        self.assertEqual(self.rescan(cve_ids=[], report=None)['action'], 'open')

    def test_a_rebuild_that_did_not_reach_quay_keeps_the_cycle_open(self):
        # The scan looks at GHCR's mega-amd64; quay.io is what users pull.
        for result in ('failure', 'skipped', 'cancelled', ''):
            with self.subTest(deploy_result=result):
                d = self.rescan(cve_ids=[], deploy_result=result, open_issue=open_issue(1150))
                self.assertEqual((d['action'], d['number']), ('comment', 1150))
                self.assertIn('quay.io', d['body'])
                self.assertIn('<!-- scan-red:deploy-failed:800 -->', d['body'])

    def test_a_rebuild_that_did_not_reach_quay_is_reported_once(self):
        noted = open_issue(1150, comments=['<!-- scan-red:deploy-failed:800 -->'])
        self.assertEqual(self.rescan(cve_ids=[], deploy_result='failure', open_issue=noted)['action'], 'none')

    def test_a_rebuild_that_did_not_reach_quay_opens_an_issue_if_none_is_open(self):
        d = self.rescan(cve_ids=[], deploy_result='failure')
        self.assertEqual(d['action'], 'open')
        self.assertIn('quay.io', d['title'])

    def test_filed_survivors_do_not_close_the_cycle_if_quay_was_not_updated(self):
        filed = [issue(1151, '[CVE-2026-7777] pillow: out-of-bounds read')]
        d = self.rescan(cve_ids=['CVE-2026-7777'], file_issues_result='success', issues=filed,
                        report=REBUILT_PILLOW, deploy_result='failure', open_issue=open_issue(1150))
        self.assertEqual(d['action'], 'comment')

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

    # --- the rescan itself fails (Trivy, the registry): an incomplete scan is never a clean one ---

    def failed_rescan(self, **kwargs):
        return self.rescan(route='none', scan_result='failure', report=None, **kwargs)

    def test_a_failed_rescan_is_reported_and_never_closes_the_cycle(self):
        for rebuild_result in ('success', 'failure'):
            with self.subTest(rebuild_result=rebuild_result):
                d = self.failed_rescan(open_issue=open_issue(1150), rebuild_result=rebuild_result)
                self.assertEqual((d['action'], d['number']), ('comment', 1150))
                for expected in (f'{REPO}/actions/runs/900', f'{REPO}/actions/runs/800', 'maintainer',
                                 '<!-- scan-red:rescan-failed:900 -->'):
                    self.assertIn(expected, d['body'])
                self.assertNotIn('cleared', d['body'])

    def test_a_failed_rescan_is_reported_once(self):
        noted = open_issue(1150, comments=['<!-- scan-red:rescan-failed:900 -->'])
        self.assertEqual(self.failed_rescan(open_issue=noted)['action'], 'none')

    def test_a_failed_rescan_opens_an_issue_if_none_is_open(self):
        d = self.failed_rescan()
        self.assertEqual(d['action'], 'open')
        self.assertIn('Rescan', d['title'])

    def test_a_failed_daily_scan_is_left_to_the_next_one(self):
        d = self.decide(route='none', scan_result='failure', open_issue=open_issue(1150))
        self.assertEqual(d['action'], 'none')

    # --- the daily scan with nothing new ---

    def test_a_daily_scan_with_nothing_new_closes_a_leftover_issue(self):
        d = self.decide(route='none', open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('No new fixable', d['body'])

    def test_a_daily_scan_with_nothing_new_and_nothing_open_does_nothing(self):
        self.assertEqual(self.decide(route='none')['action'], 'none')

    def test_a_daily_scan_with_nothing_new_leaves_open_an_issue_that_needs_a_maintainer(self):
        # A failed rebuild or a missed quay.io deploy is still there the next day.
        for marker in ('rebuild-failed:800', 'deploy-failed:800'):
            with self.subTest(marker=marker):
                stuck = open_issue(1150, body='<!-- scan-red:ids=CVE-2026-1001 -->',
                                   comments=[f'<!-- scan-red:{marker} -->'])
                d = self.decide(route='none', open_issue=stuck)
                self.assertEqual(d['action'], 'none')
                self.assertIn('maintainer', d['summary'])

    def test_a_failure_from_an_earlier_cycle_does_not_hold_a_later_one_open(self):
        later = open_issue(1150, body='<!-- scan-red:ids=CVE-2026-1001 -->',
                           comments=['<!-- scan-red:rebuild-failed:800 -->', '<!-- scan-red:ids=CVE-2026-1002 -->'])
        self.assertEqual(self.decide(route='none', open_issue=later)['action'], 'close')

    def test_nothing_new_but_a_retry_due_leaves_the_cycle_open(self):
        # e.g. a skip_rebuild run while Ubuntu's fix for an open issue is waiting for a rebuild.
        d = self.decide(route='none', retry_ids=['CVE-2026-109518'], open_issue=open_issue(1150))
        self.assertEqual(d['action'], 'none')

    def test_a_dispatched_id_still_in_the_image_below_high_is_not_called_cleared(self):
        # CVE-2026-5555 is still in the rebuilt image, as a MEDIUM finding.
        d = self.rescan(cve_ids=[], pre_rebuild_ids=['CVE-2026-109518', 'CVE-2026-5555'],
                        open_issue=open_issue(1150))
        self.assertEqual((d['action'], d['number']), ('close', 1150))
        self.assertIn('Cleared by the rebuild: CVE-2026-109518.', d['body'])
        self.assertNotIn('cleared every advisory', d['body'])
        self.assertIn('CVE-2026-5555', d['body'])
        self.assertIn('no longer a fixable HIGH/CRITICAL', d['body'])

    def test_a_direct_triage_leaves_the_issue_alone(self):
        # skip_rebuild: a maintainer sent these straight to Claude; no rebuild cycle exists.
        d = self.decide(route='triage', cve_ids=['CVE-2026-1001'], open_issue=open_issue(1150))
        self.assertEqual(d['action'], 'none')

    # --- a rebuild for an open issue whose Ubuntu fix was published ---

    def test_a_retry_rebuild_opens_an_issue_about_the_published_fix(self):
        d = self.decide(route='rebuild', retry_ids=['CVE-2026-109518'], dispatch='dispatched')
        self.assertEqual(d['action'], 'open')
        self.assertIn('Ubuntu', d['title'])
        self.assertIn('2026-11-05', d['title'])
        self.assertRegex(d['body'], r'\| CVE-2026-109518 \| HIGH \| apt, libapt-pkg7\.0 \| ubuntu \|')
        self.assertNotIn('0 new fixable', d['body'])
        self.assertIn('<!-- scan-red:ids=CVE-2026-109518 -->', d['body'])

    def test_new_and_retried_ids_share_one_cycle(self):
        d = self.decide(route='rebuild', cve_ids=['CVE-2026-7777'], retry_ids=['CVE-2026-109518'],
                        dispatch='dispatched')
        self.assertEqual(d['action'], 'open')
        self.assertIn('New fixable advisories', d['title'])
        self.assertRegex(d['body'], r'\| CVE-2026-7777 \| CRITICAL \|')
        self.assertRegex(d['body'], r'\| CVE-2026-109518 \| HIGH \|')
        self.assertIn('<!-- scan-red:ids=CVE-2026-109518,CVE-2026-7777 -->', d['body'])


BOT = {'login': 'app/github-actions', 'is_bot': True}


def cve_issue(number, title, author=BOT, labels=('security', 'cve'), comments=(), closed_at=None):
    """A `cve` issue as `gh issue list --json number,title,author,labels,comments,closedAt` gives it.
    Each comment is (login, body) or (login, body, createdAt)."""
    issue = {'number': number, 'title': title, 'author': author, 'labels': [{'name': n} for n in labels],
             'comments': [{'author': {'login': c[0]}, 'body': c[1], **({'createdAt': c[2]} if len(c) > 2 else {})}
                          for c in comments]}
    if closed_at:
        issue['closedAt'] = closed_at
    return issue


APT_ISSUE = cve_issue(1144, '[CVE-2026-109518] apt: unverified package installed as root via hard-link reuse')
PILLOW_ISSUE = cve_issue(1151, '[CVE-2026-7777] pillow: out-of-bounds read')
BROTLI_ISSUE = cve_issue(1142, '[CVE-2026-111497] libbrotli1: Out-of-bounds write in Brotli decoder')
# The apt probe's answer: {(package, fixed version): archive candidate at least that version}.
APT_FIX = {('apt', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.1', ('libapt-pkg7.0', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.1'}
# The versions request-rebuild dispatched a rebuild for, given APT_FIX.
APT_FIX_VERSIONS = ['apt=3.2.1ubuntu0.1', 'libapt-pkg7.0=3.2.1ubuntu0.1']
# The comment update-cve-issues writes when it closes an issue as gone, and the close a few seconds later.
GONE = ('github-actions', 'No longer found.\n\n<!-- cve-backlog:gone -->', '2026-11-06T09:10:02Z')
GONE_CLOSED_AT = '2026-11-06T09:10:04Z'


class BacklogTest(unittest.TestCase):
    """What the scan says about the open `cve` issues this workflow filed: gone (close it),
    waiting on Ubuntu (keep it from the fix-PR agent) and, once Ubuntu publishes the fix,
    due a rebuild."""

    def numbers(self, entries):
        return [e['number'] for e in entries]

    def test_an_issue_the_scan_no_longer_reports_is_gone(self):
        b = cp.backlog([BROTLI_ISSUE, APT_ISSUE, PILLOW_ISSUE], APT_REPORT)
        self.assertTrue(b['complete'])
        self.assertEqual(b['gone'], [{'number': 1142, 'ids': ['CVE-2026-111497']}])

    def test_an_issue_whose_findings_are_all_ubuntu_packages_is_waiting(self):
        b = cp.backlog([BROTLI_ISSUE, APT_ISSUE, PILLOW_ISSUE], APT_REPORT)
        self.assertEqual(self.numbers(b['waiting']), [1144])
        waiting = b['waiting'][0]
        self.assertEqual(waiting['vuln_ids'], ['CVE-2026-109518'])
        self.assertEqual(waiting['packages'], [{'name': 'apt', 'fixed': ['3.2.1ubuntu0.1']},
                                               {'name': 'libapt-pkg7.0', 'fixed': ['3.2.1ubuntu0.1']}])

    def test_an_issue_with_a_conda_finding_is_left_alone(self):
        b = cp.backlog([PILLOW_ISSUE], APT_REPORT)
        self.assertEqual((b['gone'], b['waiting']), ([], []))

    def test_an_issue_with_both_ubuntu_and_conda_findings_is_not_waiting(self):
        # The fix-PR agent can still fix the conda half.
        scan = {'Results': [
            {'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu',
             'Vulnerabilities': [os_finding('CVE-2026-4444', 'libssl3t64', '3.5.0-1', '3.5.0-1ubuntu0.1')]},
            {'Target': 'conda', 'Class': 'lang-pkgs', 'Type': 'conda-pkg',
             'Vulnerabilities': [os_finding('CVE-2026-4444', 'openssl', '3.5.0', '3.5.1')]},
        ]}
        b = cp.backlog([cve_issue(1160, '[CVE-2026-4444] openssl: bug')], scan)
        self.assertEqual((b['gone'], b['waiting']), ([], []))

    def test_any_finding_of_the_id_keeps_an_issue_open(self):
        # MEDIUM and unfixed findings are in trivy-results.json too; the issue is not gone.
        scan = {'Results': [{'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu'},
                            {'Target': 'jar', 'Class': 'lang-pkgs', 'Type': 'jar', 'Vulnerabilities': [
                                finding('CVE-2026-6001', severity='MEDIUM'), finding('CVE-2026-6002', fixed='')]}]}
        issues = [cve_issue(1170, '[CVE-2026-6001] foo: bar'), cve_issue(1171, '[CVE-2026-6002] foo: baz')]
        self.assertEqual(cp.backlog(issues, scan)['gone'], [])

    def test_an_issue_titled_with_an_alias_is_still_reported(self):
        # Filed under the GHSA; Trivy now reports the CVE, with the GHSA in VendorIDs.
        scan = {'Results': [{'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu'},
                            {'Target': 'py', 'Class': 'lang-pkgs', 'Type': 'python-pkg', 'Vulnerabilities': [
                                finding('CVE-2026-2002', vendor_ids=['GHSA-6q48-6c7w-9gw2'])]}]}
        ghsa_issue = cve_issue(1136, ISSUES_1137[4]['title'])
        self.assertEqual(cp.backlog([ghsa_issue], scan)['gone'], [])

    def test_only_issues_this_workflow_filed_are_managed(self):
        human = cve_issue(1180, '[CVE-2026-111497] libbrotli1: filed by hand', author={'login': 'dpark01'})
        test = cve_issue(1181, '[CVE-2026-111497] libbrotli1: test issue', labels=('security', 'cve', 'test'))
        b = cp.backlog([human, test], APT_REPORT)
        self.assertEqual((b['gone'], b['waiting']), ([], []))

    def test_an_issue_with_no_advisory_id_in_its_title_is_left_alone(self):
        b = cp.backlog([cve_issue(1182, 'Container scan follow-up')], APT_REPORT)
        self.assertEqual(b['gone'], [])

    def test_a_report_that_did_not_scan_the_os_closes_nothing(self):
        for scan in (REBUILT_PILLOW, {'Results': []}, {}):
            with self.subTest(scan=scan):
                b = cp.backlog([BROTLI_ISSUE], scan)
                self.assertFalse(b['complete'])
                self.assertEqual(b['gone'], [])

    def test_every_fixed_version_is_listed(self):
        scan = {'Results': [{'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu', 'Vulnerabilities': [
            os_finding('CVE-2026-4445', 'libfoo1', '1.0-1', '1.0-1ubuntu0.1, 1.1-1')]}]}
        b = cp.backlog([cve_issue(1190, '[CVE-2026-4445] libfoo1: bug')], scan)
        self.assertEqual(b['waiting'][0]['packages'], [{'name': 'libfoo1', 'fixed': ['1.0-1ubuntu0.1', '1.1-1']}])

    def test_names_and_versions_that_are_not_debian_shaped_are_dropped(self):
        # They reach a shell in the apt probe.
        scan = {'Results': [{'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu', 'Vulnerabilities': [
            os_finding('CVE-2026-4446', 'evil;touch x', '1.0', '1.1'),
            os_finding('CVE-2026-4446', 'libok1', '1.0', '$(touch x)'),
            os_finding('CVE-2026-4446', 'libok3', '1.0', '1.1;touch x'),
            os_finding('CVE-2026-4446', 'libok2', '1.0', '1:1.1~rc1+dfsg-2ubuntu0.1'),
        ]}]}
        b = cp.backlog([cve_issue(1191, '[CVE-2026-4446] libok: bug')], scan)
        self.assertEqual(b['waiting'][0]['packages'], [{'name': 'libok2', 'fixed': ['1:1.1~rc1+dfsg-2ubuntu0.1']}])

    # --- retrying the rebuild once Ubuntu publishes the fix ---

    def test_without_the_apt_probe_nothing_is_retried(self):
        b = cp.backlog([APT_ISSUE], APT_REPORT)
        self.assertEqual((b['retry'], b['retry_ids']), ([], []))

    def test_a_published_fix_is_due_a_rebuild(self):
        b = cp.backlog([APT_ISSUE, PILLOW_ISSUE], APT_REPORT, published=APT_FIX)
        self.assertEqual(self.numbers(b['retry']), [1144])
        self.assertEqual(b['retry'][0]['markers'], ['apt=3.2.1ubuntu0.1', 'libapt-pkg7.0=3.2.1ubuntu0.1'])
        self.assertEqual(b['retry_ids'], ['CVE-2026-109518'])

    def test_a_fix_a_rebuild_already_missed_is_not_retried(self):
        for login in ('github-actions', 'app/github-actions'):
            with self.subTest(login=login):
                tried = cve_issue(1144, APT_ISSUE['title'], comments=[
                    (login, 'Still here.\n\n<!-- cve-retry:apt=3.2.1ubuntu0.1 -->\n'
                            '<!-- cve-retry:libapt-pkg7.0=3.2.1ubuntu0.1 -->')])
                b = cp.backlog([tried], APT_REPORT, published=APT_FIX)
                self.assertEqual(self.numbers(b['waiting']), [1144])
                self.assertEqual((b['retry'], b['retry_ids']), ([], []))

    def test_a_retry_marker_from_anyone_else_is_ignored(self):
        # Anyone can comment on a public repo.
        spoofed = cve_issue(1144, APT_ISSUE['title'], comments=[
            ('mallory', '<!-- cve-retry:apt=3.2.1ubuntu0.1 --><!-- cve-retry:libapt-pkg7.0=3.2.1ubuntu0.1 -->')])
        b = cp.backlog([spoofed], APT_REPORT, published=APT_FIX)
        self.assertEqual(self.numbers(b['retry']), [1144])

    def test_a_newer_published_version_is_retried(self):
        tried = cve_issue(1144, APT_ISSUE['title'], comments=[
            ('github-actions', '<!-- cve-retry:apt=3.2.1ubuntu0.1 -->')])
        b = cp.backlog([tried], APT_REPORT, published={('apt', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.2'})
        self.assertEqual(b['retry'][0]['markers'], ['apt=3.2.1ubuntu0.2'])

    def test_a_published_version_counts_only_for_the_issues_it_fixes(self):
        # Two issues in one package: Ubuntu has published the first one's fix, not the second's.
        scan = {'Results': [{'Target': 'x (ubuntu 26.04)', 'Class': 'os-pkgs', 'Type': 'ubuntu', 'Vulnerabilities': [
            os_finding('CVE-2026-4447', 'libssl3t64', '3.5.0-1', '3.5.0-1ubuntu0.1'),
            os_finding('CVE-2026-4448', 'libssl3t64', '3.5.0-1', '3.5.0-1ubuntu0.2')]}]}
        issues = [cve_issue(1192, '[CVE-2026-4447] libssl3t64: bug'), cve_issue(1193, '[CVE-2026-4448] libssl3t64: bug')]
        b = cp.backlog(issues, scan, published={('libssl3t64', '3.5.0-1ubuntu0.1'): '3.5.0-1ubuntu0.1'})
        self.assertEqual(self.numbers(b['retry']), [1192])
        self.assertEqual(b['retry_ids'], ['CVE-2026-4447'])

    # --- reopening what the scan closed, when the advisory comes back ---

    def closed(self, title='[CVE-2026-7777] pillow: out-of-bounds read', comments=(GONE,), closed_at=GONE_CLOSED_AT,
               **kwargs):
        return cve_issue(1151, title, comments=comments, closed_at=closed_at, **kwargs)

    def reopened(self, issue):
        return cp.backlog([], APT_REPORT, closed_issues=[issue])['reopen']

    def test_an_issue_the_scan_closed_is_reopened_when_its_advisory_returns(self):
        # e.g. a .trivyignore.yaml entry expired, or a Trivy DB gap closed it for a day.
        b = cp.backlog([], APT_REPORT, closed_issues=[self.closed()])
        self.assertEqual(b['reopen'], [{'number': 1151, 'ids': ['CVE-2026-7777']}])

    def test_an_issue_closed_any_other_way_stays_closed(self):
        # A fix PR or a maintainer closed it: the scan does not second-guess that.
        for comments in ((), [('dpark01', 'Fixed by #1137', '2026-11-06T09:10:02Z')]):
            with self.subTest(comments=comments):
                self.assertEqual(self.reopened(self.closed(comments=comments)), [])

    def test_an_issue_closed_again_after_a_reopen_stays_closed(self):
        reclosed = self.closed(comments=[GONE, ('github-actions', '<!-- cve-backlog:reopened -->', '2026-11-07T09:10:02Z'),
                                         ('dpark01', 'Accepted risk, closing.', '2026-11-07T10:00:00Z')],
                               closed_at='2026-11-07T10:00:01Z')
        self.assertEqual(self.reopened(reclosed), [])

    def test_a_maintainer_who_reopens_and_closes_again_is_not_overruled(self):
        # Reopened and closed again by hand, with no comment at all: the latest close is the
        # maintainer's, not the scan's, though the scan's gone comment is still the latest marker.
        self.assertEqual(self.reopened(self.closed(closed_at='2026-11-08T15:00:00Z')), [])

    def test_a_close_before_the_gone_comment_is_not_the_scans(self):
        self.assertEqual(self.reopened(self.closed(closed_at='2026-11-06T09:00:00Z')), [])

    def test_without_the_close_and_comment_times_nothing_is_reopened(self):
        for issue in (self.closed(closed_at=None),
                      self.closed(comments=[('github-actions', '<!-- cve-backlog:gone -->')]),
                      self.closed(closed_at='yesterday')):
            with self.subTest(issue=issue):
                self.assertEqual(self.reopened(issue), [])

    def test_a_gone_marker_from_anyone_else_is_ignored(self):
        spoofed = self.closed(comments=[('mallory', '<!-- cve-backlog:gone -->', '2026-11-06T09:10:02Z')])
        self.assertEqual(self.reopened(spoofed), [])

    def test_an_advisory_still_gone_stays_closed(self):
        self.assertEqual(self.reopened(self.closed(title='[CVE-2026-111497] libbrotli1: bug')), [])

    def test_only_issues_this_workflow_filed_are_reopened(self):
        for kwargs in ({'author': {'login': 'dpark01'}}, {'labels': ('security', 'cve', 'test')}):
            with self.subTest(**kwargs):
                self.assertEqual(self.reopened(self.closed(**kwargs)), [])

    # --- what the rescan and update-cve-issues need from the backlog ---

    def test_the_backlog_records_the_digest_it_scanned(self):
        self.assertEqual(cp.backlog([], APT_REPORT)['scanned_digest'], MAIN_DIGEST)
        for metadata in ({}, {'RepoDigests': []}, {'RepoDigests': ['ghcr.io/broadinstitute/viral-ngs@sha256:xyz']}):
            with self.subTest(metadata=metadata):
                self.assertIsNone(cp.backlog([], {**APT_REPORT, 'Metadata': metadata})['scanned_digest'])

    def test_the_backlog_lists_the_versions_a_rebuild_is_dispatched_for(self):
        b = cp.backlog([APT_ISSUE], APT_REPORT, published=APT_FIX)
        self.assertEqual(b['retry_versions'], APT_FIX_VERSIONS)
        self.assertEqual(cp.backlog([APT_ISSUE], APT_REPORT)['retry_versions'], [])

    def test_only_waiting_issues_are_retried(self):
        b = cp.backlog([PILLOW_ISSUE], APT_REPORT, published={('pillow', '13.0.0'): '13.0.0'})
        self.assertEqual(b['retry'], [])


class BacklogActionsTest(unittest.TestCase):

    def actions(self, backlog, rebuild_state='not-rebuilt', post_rebuild_run_id='', pre_rebuild_versions=(),
                main_digest=MAIN_DIGEST):
        return cp.backlog_actions(backlog, rebuild_state=rebuild_state, post_rebuild_run_id=post_rebuild_run_id,
                                  pre_rebuild_versions=pre_rebuild_versions, main_digest=main_digest, repo_url=REPO,
                                  run_id='900', image='ghcr.io/broadinstitute/viral-ngs@sha256:abc')

    def rebuilt(self, backlog, pre_rebuild_versions=APT_FIX_VERSIONS, **kwargs):
        return self.actions(backlog, rebuild_state='rebuilt', post_rebuild_run_id='800',
                            pre_rebuild_versions=pre_rebuild_versions, **kwargs)

    def test_a_gone_issue_is_closed_with_the_run_that_saw_it_gone(self):
        [action] = self.actions(cp.backlog([BROTLI_ISSUE, PILLOW_ISSUE], APT_REPORT))
        self.assertEqual((action['action'], action['number']), ('close', 1142))
        for expected in ('CVE-2026-111497', f'{REPO}/actions/runs/900', 'sha256:abc', '<!-- cve-backlog:gone -->'):
            self.assertIn(expected, action['body'])

    def test_a_returning_advisory_reopens_its_issue(self):
        gone = cve_issue(1151, PILLOW_ISSUE['title'], comments=[GONE], closed_at=GONE_CLOSED_AT)
        [action] = self.actions(cp.backlog([], APT_REPORT, closed_issues=[gone]))
        self.assertEqual((action['action'], action['number']), ('reopen', 1151))
        for expected in ('CVE-2026-7777', f'{REPO}/actions/runs/900', '<!-- cve-backlog:reopened -->'):
            self.assertIn(expected, action['body'])

    def test_nothing_is_closed_or_reopened_once_main_has_moved_on(self):
        # A push since the scan (or since the rebuild it rescans) published a newer image the
        # scan never saw; an unresolvable tag counts as moved. The next daily scan decides.
        gone = cve_issue(1151, PILLOW_ISSUE['title'], comments=[GONE], closed_at=GONE_CLOSED_AT)
        b = cp.backlog([BROTLI_ISSUE], APT_REPORT, closed_issues=[gone])
        self.assertEqual(len(self.actions(b)), 2)
        for main_digest in ('sha256:' + 'b' * 64, ''):
            with self.subTest(main_digest=main_digest):
                self.assertEqual(self.actions(b, main_digest=main_digest), [])
        unpinned = cp.backlog([BROTLI_ISSUE], {**APT_REPORT, 'Metadata': {}})
        self.assertEqual(self.actions(unpinned), [])

    def test_a_missed_fix_is_still_marked_once_main_has_moved_on(self):
        # The rebuild's image did not get the fix whatever main points at now.
        b = cp.backlog([APT_ISSUE], APT_REPORT, published=APT_FIX)
        [action] = self.rebuilt(b, main_digest='sha256:' + 'b' * 64)
        self.assertEqual((action['action'], action['number']), ('comment', 1144))

    def test_a_rescan_marks_a_published_fix_its_rebuild_did_not_install(self):
        b = cp.backlog([APT_ISSUE], APT_REPORT, published=APT_FIX)
        [action] = self.rebuilt(b)
        self.assertEqual((action['action'], action['number']), ('comment', 1144))
        for expected in (f'{REPO}/actions/runs/800', 'maintainer', '3.2.1ubuntu0.1',
                         '<!-- cve-retry:apt=3.2.1ubuntu0.1 -->', '<!-- cve-retry:libapt-pkg7.0=3.2.1ubuntu0.1 -->'):
            self.assertIn(expected, action['body'])

    def test_only_a_successful_rebuild_marks_a_fix_as_tried(self):
        # Before a rebuild (the route dispatches one) or after a failed one, the fix is still untried.
        b = cp.backlog([APT_ISSUE], APT_REPORT, published=APT_FIX)
        for state in ('not-rebuilt', 'rebuild-failed'):
            with self.subTest(rebuild_state=state):
                self.assertEqual(self.actions(b, rebuild_state=state, post_rebuild_run_id='800',
                                              pre_rebuild_versions=APT_FIX_VERSIONS), [])

    def test_only_the_fixes_the_rebuild_was_dispatched_for_are_marked(self):
        # The rebuild ran for another package's fix, or was not the scan's (no versions passed).
        b = cp.backlog([APT_ISSUE], APT_REPORT, published=APT_FIX)
        for versions in (['libpcre2-8-0=10.46-1ubuntu0.1'], []):
            with self.subTest(pre_rebuild_versions=versions):
                self.assertEqual(self.rebuilt(b, pre_rebuild_versions=versions), [])

    def test_a_version_published_after_the_dispatch_is_left_for_another_rebuild(self):
        # The rebuild was dispatched for 0.1; Ubuntu published 0.2 before the rescan's probe ran.
        # Only what the rebuild could have installed counts as tried.
        newer = {('apt', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.2', ('libapt-pkg7.0', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.2'}
        self.assertEqual(self.rebuilt(cp.backlog([APT_ISSUE], APT_REPORT, published=newer)), [])

        partly = {('apt', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.1', ('libapt-pkg7.0', '3.2.1ubuntu0.1'): '3.2.1ubuntu0.2'}
        [action] = self.rebuilt(cp.backlog([APT_ISSUE], APT_REPORT, published=partly))
        self.assertIn('<!-- cve-retry:apt=3.2.1ubuntu0.1 -->', action['body'])
        self.assertNotIn('cve-retry:libapt-pkg7.0', action['body'])
        self.assertNotIn('3.2.1ubuntu0.2', action['body'])

    def test_nothing_to_do(self):
        self.assertEqual(self.actions(cp.backlog([PILLOW_ISSUE], APT_REPORT)), [])


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
            (['--cve-ids', '', '--retry-ids', 'CVE-2026-109518'], 'rebuild'),
            (['--cve-ids', '', '--retry-ids', 'CVE-2026-109518', '--post-rebuild-run-id', '800'], 'none'),
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
            '--deploy-result', 'success',
            '--dispatch', '',
            '--file-issues-result', 'skipped',
            '--open-issue-json', self.path('open.json', open_issue(1150)),
            '--issues-json', self.path('issues.json', []),
            '--trivy-json', self.path('trivy.json', REBUILT_CLEAN),
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
            '--issues-json', self.path('issues.json', []),
            '--trivy-json', self.path('missing-trivy.json'),
            '--repo-url', REPO,
            '--run-id', '900',
            '--image', 'ghcr.io/broadinstitute/viral-ngs:main-mega-amd64',
            '--date', '2026-11-05',
            '--out', self.path('decision.json'),
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(self.read('decision.json'))['action'], 'open')

    def test_scan_red_for_a_retry(self):
        rc = self.run_main([
            'scan-red',
            '--route', 'rebuild',
            '--cve-ids', '',
            '--retry-ids', 'CVE-2026-109518',
            '--dispatch', 'dispatched',
            '--open-issue-json', self.path('open.json', 'null'),
            '--issues-json', self.path('issues.json', []),
            '--trivy-json', self.path('trivy.json', APT_REPORT),
            '--repo-url', REPO,
            '--run-id', '900',
            '--image', 'ghcr.io/broadinstitute/viral-ngs:main-mega-amd64',
            '--date', '2026-11-05',
            '--out', self.path('decision.json'),
        ])
        self.assertEqual(rc, 0)
        decision = json.loads(self.read('decision.json'))
        self.assertEqual(decision['action'], 'open')
        self.assertIn('<!-- scan-red:ids=CVE-2026-109518 -->', decision['body'])

    def test_backlog_lists_the_packages_to_probe_then_the_retries(self):
        issues = self.path('issues.json', [BROTLI_ISSUE, APT_ISSUE, PILLOW_ISSUE])
        trivy = self.path('trivy.json', APT_REPORT)
        rc = self.run_main(['backlog', '--trivy-json', trivy, '--issues-json', issues,
                            '--waiting-pkgs-out', self.path('waiting.txt'), '--out', self.path('backlog.json')])
        self.assertEqual(rc, 0)
        self.assertEqual(self.read('waiting.txt'), 'apt 3.2.1ubuntu0.1\nlibapt-pkg7.0 3.2.1ubuntu0.1\n')
        self.assertEqual(json.loads(self.read('backlog.json'))['retry_ids'], [])

        # The probe prints "<package> <fixed> <candidate>" for each fixed version that is published;
        # a malformed line is dropped.
        published = self.path('published.txt', 'apt 3.2.1ubuntu0.1 3.2.1ubuntu0.1\n'
                              'libapt-pkg7.0 3.2.1ubuntu0.1 3.2.1ubuntu0.1\nbad line\n')
        rc = self.run_main(['backlog', '--trivy-json', trivy, '--issues-json', issues, '--published', published,
                            '--out', self.path('backlog.json')])
        self.assertEqual(rc, 0)
        b = json.loads(self.read('backlog.json'))
        self.assertEqual(b['retry_ids'], ['CVE-2026-109518'])
        self.assertEqual([g['number'] for g in b['gone']], [1142])

        gone = cve_issue(1151, PILLOW_ISSUE['title'], comments=[GONE], closed_at=GONE_CLOSED_AT)
        rc = self.run_main(['backlog', '--trivy-json', trivy, '--issues-json', issues,
                            '--closed-issues-json', self.path('closed.json', [gone]), '--out', self.path('backlog.json')])
        self.assertEqual(rc, 0)
        b = json.loads(self.read('backlog.json'))
        self.assertEqual([r['number'] for r in b['reopen']], [1151])
        self.assertEqual(b['scanned_digest'], MAIN_DIGEST)

    def backlog_actions(self, *extra):
        backlog = self.path('backlog.json', cp.backlog([BROTLI_ISSUE, APT_ISSUE], APT_REPORT, published=APT_FIX))
        with mock.patch('sys.stdout', io.StringIO()) as out:
            rc = cp.main(['backlog-actions', '--backlog', backlog, '--rebuild-state', 'rebuilt',
                          '--post-rebuild-run-id', '800', '--repo-url', REPO, '--run-id', '900',
                          '--image', 'ghcr.io/broadinstitute/viral-ngs@sha256:abc',
                          '--out', self.path('actions.json'), *extra])
        self.assertEqual(rc, 0)
        return [(a['action'], a['number']) for a in json.loads(self.read('actions.json'))], out.getvalue()

    def test_backlog_actions(self):
        # Malformed versions are dropped, not passed into issue text.
        actions, _ = self.backlog_actions('--pre-rebuild-versions', ' '.join(APT_FIX_VERSIONS) + ' x;y=1 apt=1;id',
                                          '--main-digest', MAIN_DIGEST)
        self.assertEqual(actions, [('close', 1142), ('comment', 1144)])
        actions, _ = self.backlog_actions('--pre-rebuild-versions', 'apt=1;id libapt-pkg7.0=3.2.1ubuntu0.1;id',
                                          '--main-digest', MAIN_DIGEST)
        self.assertEqual(actions, [('close', 1142)])

    def test_backlog_actions_once_main_has_moved_on(self):
        actions, out = self.backlog_actions('--pre-rebuild-versions', ' '.join(APT_FIX_VERSIONS),
                                            '--main-digest', '')
        self.assertEqual(actions, [('comment', 1144)])
        self.assertIn('::warning::', out)
        self.assertIn('#1142', out)

    def test_scan_red_for_a_failed_rescan(self):
        rc = self.run_main([
            'scan-red', '--route', 'none', '--scan-result', 'failure', '--cve-ids', '',
            '--post-rebuild-run-id', '800', '--rebuild-result', 'success', '--deploy-result', 'success',
            '--open-issue-json', self.path('open.json', open_issue(1150)),
            '--issues-json', self.path('issues.json', []),
            '--trivy-json', self.path('missing.json'),
            '--repo-url', REPO, '--run-id', '900', '--image', 'ghcr.io/broadinstitute/viral-ngs@sha256:abc',
            '--date', '2026-11-06', '--out', self.path('decision.json'),
        ])
        self.assertEqual(rc, 0)
        decision = json.loads(self.read('decision.json'))
        self.assertEqual((decision['action'], decision['number']), ('comment', 1150))


if __name__ == '__main__':
    unittest.main()
