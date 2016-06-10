# Unit tests for metagenomics.py
from builtins import super
import six
import argparse
from collections import Counter
import copy
import io
import os.path
from os.path import join
import tempfile
import textwrap
import unittest
import pytest

import mock
from mock import patch

import tools.picard
import metagenomics
import util.file
import util.misc
from test import TestCaseWithTmp


if six.PY2:
    from StringIO import StringIO
else:
    from io import StringIO


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in metagenomics.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestKrakenCalls(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = patch('tools.picard.SamToFastqTool', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_samtofastq = patcher.start()

        patcher = patch('tools.kraken.Kraken', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_kraken = patcher.start()

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_assert(self):
        with self.assertRaises(AssertionError):
            metagenomics.kraken(self.inBam, self.db)

    def test_out_reads(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads)
        self.mock_samtofastq().execute.assert_called_once_with(
            self.inBam, mock.ANY, mock.ANY, picardOptions=mock.ANY,
            JVMmemory=mock.ANY)
        picard_opts = self.mock_samtofastq().execute.call_args[1]['picardOptions']
        six.assertCountEqual(
            self, picard_opts,
            ['CLIPPING_ACTION=X', 'CLIPPING_ATTRIBUTE=%s'
             % tools.picard.SamToFastqTool.illumina_clipping_attribute])

        self.mock_kraken().classify.assert_called_once_with(
            mock.ANY, self.db, mock.ANY, mock.ANY)

        # Exists due to test. Need a better tmpfile patcher.
        self.assertTrue(os.path.isfile(out_reads))

    def test_filter_threshold(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, filterThreshold=0.05)
        self.mock_kraken().filter.assert_called_with(
            mock.ANY, self.db, mock.ANY, 0.05)
        self.mock_kraken().classify.assert_called_once_with(
            mock.ANY, self.db, mock.ANY, mock.ANY)

    def test_out_report(self):
        out_report = util.file.mkstempfname('.kraken_report.txt')
        metagenomics.kraken(self.inBam, self.db, outReport=out_report)

        self.mock_kraken().classify.assert_called_once_with(
            mock.ANY, self.db, mock.ANY, mock.ANY)

        self.mock_kraken().report.assert_called_with(
            mock.ANY, self.db, out_report)

    def test_out_reads_and_report(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        out_report = util.file.mkstempfname('.kraken_report.txt')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, outReport=out_report)
        self.mock_kraken().classify.assert_called_once_with(
            self.inBam, self.db, mock.ANY, mock.ANY)
        self.mock_kraken().report.assert_called_with(
            mock.ANY, self.db, out_report)

    def test_num_threads(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, numThreads=11)
        self.mock_kraken().classify.assert_called_once_with(
            mock.ANY, self.db, mock.ANY, numThreads=min(util.misc.available_cpu_count(), 11)
        )


@patch('metagenomics.blast_lca')
@patch('metagenomics.kraken_dfs_report')
class TestDiamondCalls(TestCaseWithTmp):
    def setUp(self):
        super().setUp()
        patcher = patch('tools.picard.SamToFastqTool', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_samtofastq = patcher.start()

        patcher = patch('tools.diamond.Diamond', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_diamond = patcher.start()

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')
        self.tax_db = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple', 'db', 'taxonomy')

    def test_output_lca(self, mock_blast_lca, mock_dfs):
        out_report = util.file.mkstempfname('report.txt')
        metagenomics.diamond(self.inBam, self.db, self.tax_db, out_report)
        self.mock_samtofastq().execute.assert_called_once_with(
            self.inBam, mock.ANY, mock.ANY, picardOptions=mock.ANY,
            JVMmemory=mock.ANY)
        picard_opts = self.mock_samtofastq().execute.call_args[1]['picardOptions']
        six.assertCountEqual(
            self, picard_opts,
            ['CLIPPING_ACTION=X', 'CLIPPING_ATTRIBUTE=%s'
             % tools.picard.SamToFastqTool.illumina_clipping_attribute])
        self.mock_diamond().install.assert_called_once_with()
        self.mock_diamond().blastx.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={'--threads': 1})
        assert self.mock_diamond().view.called
        assert isinstance(metagenomics.blast_lca.call_args[0][0], metagenomics.TaxonomyDb)
        assert isinstance(metagenomics.kraken_dfs_report.call_args[0][0], metagenomics.TaxonomyDb)

    def test_num_threads(self, mock_blast_lca, mock_dfs):
        out_report = util.file.mkstempfname('report.txt')
        metagenomics.diamond(self.inBam, self.db, self.tax_db, out_report, numThreads=11)
        self.mock_diamond().blastx.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={'--threads': 11})
        self.assertEqual(self.mock_diamond().view.call_args[1]['options']['--threads'], 11)


class TestKronaCalls(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = patch('tools.krona.Krona', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_krona = patcher.start()

        self.inTsv = util.file.mkstempfname('.tsv')
        self.db = tempfile.mkdtemp('db')

    def test_krona_import_taxonomy(self):
        out_html = util.file.mkstempfname('.html')
        metagenomics.krona(self.inTsv, self.db, out_html, queryColumn=3, taxidColumn=5, scoreColumn=7,
                           noHits=True, noRank=True)
        self.mock_krona().import_taxonomy.assert_called_once_with(
            self.db, [self.inTsv], out_html, query_column=3, taxid_column=5, score_column=7,
            no_hits=True, no_rank=True)


@pytest.fixture
def taxa_db_simple():
    db = metagenomics.TaxonomyDb()
    db.gis = {1:2, 2:3, 3:4, 4:5}
    db.parents = {1: 1, 2: 1, 3: 2, 4: 3, 5: 4}
    return db


@pytest.fixture
def taxa_db(parents, names, ranks):
    db = metagenomics.TaxonomyDb()
    db.parents = parents
    db.names = names
    db.ranks = ranks
    return db


@pytest.fixture
def parents():
    return {
        1: 1,
        3: 1,
        6: 3,
        7: 6,
        8: 7,
        10: 6,
        11: 7,
        12: 8,
        13: 12
    }


@pytest.fixture
def names():
    names = {
        1: "root",
        2: "two",
        3: "three",
        4: "four",
        6: "six",
        7: "seven",
        8: "eight",
        10: "ten",
        11: "eleven",
        12: "twelve",
        13: "thirteen"
    }
    return names


@pytest.fixture
def ranks():
    return {
        1: "root",
        2: "",
        3: "superkingdom",
        4: "kingdom",
        6: "",
        7: "",
        8: "",
        10: "",
        11: "genus",
        12: "species",
        13: ""
    }


@pytest.fixture
def simple_m8():
    test_path = join(util.file.get_test_input_path(),
                             'TestTaxonomy')
    return open(join(test_path, 'simple.m8'))


def test_tree_level_lookup(parents):
    level_cache = {1: 1}
    assert metagenomics.tree_level_lookup(parents, 1, level_cache) == 1
    assert metagenomics.tree_level_lookup(parents, 3, level_cache) == 2
    assert metagenomics.tree_level_lookup(parents, 12, level_cache) == 6
    level_cache = {1: 1}
    assert metagenomics.tree_level_lookup(parents, 12, level_cache) == 6
    assert metagenomics.tree_level_lookup(parents, 8, level_cache) == 5


def test_push_up_tree_hits(parents):
    hits = Counter({1: 3, 3: 5, 6: 3, 7: 3, 13: 5})
    with pytest.raises(AssertionError):
        metagenomics.push_up_tree_hits(parents, hits)

    expected = hits.copy()
    assert metagenomics.push_up_tree_hits(parents, hits.copy(), min_support=1) == expected

    expected = Counter({3: 5, 6: 6, 13: 5})
    assert metagenomics.push_up_tree_hits(parents, hits.copy(), min_support=5) == expected

    assert (metagenomics.push_up_tree_hits(parents, hits.copy(), min_support=10) ==
            Counter({6: 11}))
    assert (metagenomics.push_up_tree_hits(parents, hits.copy(), min_support=18) ==
            Counter({1: 19}))
    assert (metagenomics.push_up_tree_hits(parents, hits.copy(), min_support_percent=100) ==
            Counter({1: 19}))


def test_parents_to_children(parents):
    children = metagenomics.parents_to_children(parents)
    assert children[1] == [3]


def test_rank_code():
    assert metagenomics.rank_code('species') == 'S'
    assert metagenomics.rank_code('genus') == 'G'
    assert metagenomics.rank_code('superkingdom') == 'D'
    assert metagenomics.rank_code('nonexist') == '-'


def test_blast_records(simple_m8):
    test_path = join(util.file.get_test_input_path(),
                             'TestTaxonomy')
    with simple_m8 as f:
        records = list(metagenomics.blast_records(f))
    assert len(records) == 110
    assert records[0].bit_score == 63.5
    assert records[-1].bit_score == 67.4


def test_blast_lca(taxa_db_simple, simple_m8):
    test_path = join(util.file.get_test_input_path(),
                             'TestTaxonomy')
    expected = textwrap.dedent("""\
    M04004:13:000000000-AGV3H:1:1101:12068:2105\t2
    M04004:13:000000000-AGV3H:1:1101:13451:2146\t2
    M04004:13:000000000-AGV3H:1:1101:13509:2113\t2
    M04004:13:000000000-AGV3H:1:1101:14644:2160\t2
    M04004:13:000000000-AGV3H:1:1101:18179:2130\t2
    M04004:13:000000000-AGV3H:1:1111:10629:2610\t2
    M04004:13:000000000-AGV3H:1:1111:10629:26101\t2
    """)
    out = StringIO()
    with simple_m8 as f:
        metagenomics.blast_lca(taxa_db_simple, f, out, paired=True)
        out.seek(0)
        assert out.read() == expected


def test_paired_query_id():
    tup = ['query', 'gi|10|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5]

    blast1 = metagenomics.BlastRecord(*tup)
    assert metagenomics.paired_query_id(blast1) == blast1

    new_tup = copy.copy(tup)
    new_tup[0] = 'query/1'
    new_blast1 = metagenomics.BlastRecord(*new_tup)
    assert metagenomics.paired_query_id(new_blast1) == blast1

    new_tup = copy.copy(tup)
    new_tup[0] = 'query/2'
    new_blast1 = metagenomics.BlastRecord(*new_tup)
    assert metagenomics.paired_query_id(new_blast1) == blast1

    new_tup = copy.copy(tup)
    new_tup[0] = 'query/3'
    new_blast1 = metagenomics.BlastRecord(*new_tup)
    assert metagenomics.paired_query_id(new_blast1) == new_blast1


def test_translate_gi_to_tax_id(taxa_db_simple):
    tup = ['query', 'gi|4|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5]
    blast1 = metagenomics.BlastRecord(*tup)

    tup[1] = 5
    expected = metagenomics.BlastRecord(*tup)
    assert metagenomics.translate_gi_to_tax_id(taxa_db_simple, blast1) == expected


def test_kraken_dfs_report(taxa_db):
    hits = Counter({1: 101, 3: 103, 10: 105, 12: 107})

    expected = textwrap.dedent('''\
    100.00\t416\t101\t-\t1\troot
    75.72\t315\t103\tD\t3\t  three
    50.96\t212\t0\t-\t6\t    six
    25.24\t105\t105\t-\t10\t      ten
    25.72\t107\t0\t-\t7\t      seven
    25.72\t107\t0\t-\t8\t        eight
    25.72\t107\t107\tS\t12\t          twelve
    ''')
    report = metagenomics.kraken_dfs_report(taxa_db, hits)
    text_report = '\n'.join(list(report)) + '\n'
    assert text_report == expected
