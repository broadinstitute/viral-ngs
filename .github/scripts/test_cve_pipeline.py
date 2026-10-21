"""Unit tests for cve_pipeline.py (no network access).

Run with: python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
"""

import io
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


if __name__ == '__main__':
    unittest.main()
