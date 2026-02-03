# Unit tests for krona
from builtins import super
import os.path
import tempfile
import unittest
from mock import patch
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc
from viral_ngs.classify import krona
from tests import TestCaseWithTmp


class TestToolKrona(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        self.krona = krona.Krona()
        self.krona.install()

        patcher = patch('subprocess.check_call', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run = patcher.start()

        self.db = tempfile.mkdtemp('db')

    def test_import_taxonomy(self):
        in_tsv = util_file.mkstempfname('.tsv')
        output = util_file.mkstempfname('.output')
        self.krona.import_taxonomy(
            self.db, [in_tsv], output, query_column=3, taxid_column=5, score_column=7, no_hits=True, no_rank=True
        )
        args = self.mock_run.call_args[0][0]
        self.assertEqual('ktImportTaxonomy', os.path.basename(args[0]))
        self.assertTrue(util_misc.list_contains(['-tax', self.db], args))
        self.assertTrue(util_misc.list_contains(['-q', '3'], args))
        self.assertTrue(util_misc.list_contains(['-t', '5'], args))
        self.assertTrue(util_misc.list_contains(['-s', '7'], args))
        self.assertIn('-i', args)
        self.assertIn('-k', args)

        self.krona.import_taxonomy(self.db, [in_tsv], output)
        args = self.mock_run.call_args[0][0]
        self.assertEqual('ktImportTaxonomy', os.path.basename(args[0]))
        self.assertTrue(util_misc.list_contains(['-tax', self.db], args))
        self.assertFalse(util_misc.list_contains(['-q', '3'], args))
        self.assertFalse(util_misc.list_contains(['-t', '5'], args))
        self.assertFalse(util_misc.list_contains(['-s', '7'], args))
        self.assertNotIn('-i', args)
        self.assertNotIn('-k', args)

    def test_create_db(self):
        self.krona.create_db(self.db)
        args = self.mock_run.call_args[0][0]
        self.assertEqual('updateTaxonomy.sh', os.path.basename(args[0]))
        self.assertTrue(util_misc.list_contains(['--only-build', self.db], args))
        self.assertIn(self.db, args)
