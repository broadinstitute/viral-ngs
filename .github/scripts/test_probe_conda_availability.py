"""Unit tests for probe-conda-availability.py's Trivy lookup (no network access).

Run with: python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
"""

import importlib.util
import json
import os
import tempfile
import unittest

_spec = importlib.util.spec_from_file_location(
    'probe_conda_availability', os.path.join(os.path.dirname(__file__), 'probe-conda-availability.py'))
probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe)


class CollectTargetsTest(unittest.TestCase):

    def collect(self, vulns, ids):
        with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
            json.dump({'Results': [{'Target': 'image', 'Class': 'lang-pkgs', 'Type': 'python-pkg',
                                    'Vulnerabilities': vulns}]}, f)
        self.addCleanup(os.unlink, f.name)
        errors = []
        return probe.collect_targets(f.name, ids, errors), errors

    def test_finding_is_found_by_its_vulnerability_id(self):
        targets, errors = self.collect(
            [{'VulnerabilityID': 'CVE-2026-112403', 'PkgName': 'pillow', 'FixedVersion': '13.0.0'}],
            ['CVE-2026-112403'])
        self.assertEqual(list(targets), ['pillow'])
        self.assertEqual(targets['pillow']['cve_ids'], ['CVE-2026-112403'])
        self.assertEqual(targets['pillow']['fixed_versions'], ['13.0.0'])
        self.assertEqual(errors, [])

    def test_issue_filed_under_a_ghsa_finds_the_finding_after_trivy_rekeys_it_to_a_cve(self):
        # The auto-fix job takes IDs from issue titles. An issue filed under a GHSA ID keeps
        # it after Trivy starts reporting the finding as a CVE, with the GHSA in VendorIDs.
        targets, errors = self.collect(
            [{'VulnerabilityID': 'CVE-2026-999999', 'VendorIDs': ['GHSA-6q48-6c7w-9gw2'],
              'PkgName': 'certifi', 'FixedVersion': '2026.10.19'}],
            ['GHSA-6q48-6c7w-9gw2'])
        self.assertEqual(list(targets), ['certifi'])
        self.assertEqual(targets['certifi']['cve_ids'], ['GHSA-6Q48-6C7W-9GW2'])
        self.assertEqual(targets['certifi']['reported_as'], ['CVE-2026-999999'])
        self.assertEqual(errors, [])

    def test_requested_id_absent_from_the_scan_is_reported(self):
        targets, errors = self.collect(
            [{'VulnerabilityID': 'CVE-2026-112403', 'PkgName': 'pillow', 'FixedVersion': '13.0.0'}],
            ['CVE-2026-112404'])
        self.assertEqual(targets, {})
        self.assertEqual(len(errors), 1)
        self.assertIn('CVE-2026-112404', errors[0])


if __name__ == '__main__':
    unittest.main()
