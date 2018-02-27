# Unit tests for metagenomics.py
from builtins import super
import six
import argparse
from collections import Counter
import copy
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
from test import TestCaseWithTmp, _CPUS

if six.PY2:
    from StringIO import StringIO
else:
    from io import StringIO


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in metagenomics.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


@patch('metagenomics.kraken_dfs_report')
class TestDiamondCalls(TestCaseWithTmp):
    def setUp(self):
        super().setUp()
        patcher = patch('subprocess.Popen')
        self.addCleanup(patcher.stop)
        self.mock_popen = patcher.start()
        self.mock_popen.return_value.returncode = 0

        patcher = patch('tools.picard.SamToFastqTool')
        self.addCleanup(patcher.stop)
        self.mock_s2f = patcher.start()
        self.mock_s2f.return_value.execute.return_value.returncode = 0

        patcher = patch('tools.diamond.Diamond', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_diamond = patcher.start()

        # Can't open unwritten named pipes
        if six.PY2:
            patcher = patch('__builtin__.open', mock.mock_open(read_data="id1\t1\n"))
        else:
            patcher = patch('builtins.open', mock.mock_open(read_data="id1\t1\n"))
        self.addCleanup(patcher.stop)
        patcher.start()

        # mock_open doesn't have __next__ for csv.reader
        patcher = patch('metagenomics.taxa_hits_from_tsv', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_taxa_hits = patcher.start()
        self.mock_taxa_hits.return_value = Counter({1: 100, 2: 40})

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')
        self.tax_db = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple', 'db', 'taxonomy')

    def test_output_reads(self, mock_dfs):
        out_report = util.file.mkstempfname('report.txt')
        out_reads = util.file.mkstempfname('lca.gz')

        metagenomics.diamond(self.inBam, self.db, self.tax_db, out_report, outReads=out_reads)

        cmd = self.mock_popen.call_args[0][0]
        self.assertIn(out_reads, cmd)
        assert isinstance(metagenomics.kraken_dfs_report.call_args[0][0], metagenomics.TaxonomyDb)

    def test_num_threads(self, mock_dfs):
        out_report = util.file.mkstempfname('report.txt')
        metagenomics.diamond(self.inBam, self.db, self.tax_db, out_report, threads=11)
        expected_threads = min(11, _CPUS)
        expected_threads = '--threads {}'.format(expected_threads)
        cmd = self.mock_popen.call_args[0][0]
        self.assertIn(expected_threads, cmd)


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
            no_hits=True, no_rank=True, magnitude_column=None, root_name=os.path.basename(self.inTsv))


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
    expected = [
        ('M04004:13:000000000-AGV3H:1:1101:12068:2105', 2),
        ('M04004:13:000000000-AGV3H:1:1101:13451:2146', 2),
        ('M04004:13:000000000-AGV3H:1:1101:13509:2113', 2),
        ('M04004:13:000000000-AGV3H:1:1101:14644:2160', 2),
        ('M04004:13:000000000-AGV3H:1:1101:18179:2130', 2),
        ('M04004:13:000000000-AGV3H:1:1111:10629:2610', 2),
        ('M04004:13:000000000-AGV3H:1:1111:10629:26101', 2),
    ]
    with simple_m8 as f:
        out = StringIO()
        hits = list(metagenomics.blast_lca(taxa_db_simple, f, out, paired=True))
        assert hits == expected


def test_paired_query_id():
    tup = ['query', 'gi|10|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5, []]

    blast1 = metagenomics.BlastRecord(*tup)
    assert metagenomics.paired_query_id(blast1) == blast1

    new_blast1 = copy.copy(blast1)
    new_blast1.query_id = 'query/1'
    assert metagenomics.paired_query_id(new_blast1) == blast1

    new_blast1 = copy.copy(blast1)
    new_blast1.query_id = 'query/2'
    assert metagenomics.paired_query_id(new_blast1) == blast1

    new_blast1 = copy.copy(blast1)
    new_blast1.query_id = 'query/3'
    assert metagenomics.paired_query_id(new_blast1) == new_blast1


def test_fill_tax_id_from_gi(taxa_db_simple):
    tup = ['query', 'gi|4|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5, []]
    blast1 = metagenomics.BlastRecord(*tup)

    metagenomics.fill_tax_id_from_gi(taxa_db_simple, blast1)
    assert blast1.taxids == [5]


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


def test_coverage_lca(taxa_db):
    assert metagenomics.coverage_lca([10, 11, 12], taxa_db.parents) == 6
    assert metagenomics.coverage_lca([1, 3], taxa_db.parents) == 1
    assert metagenomics.coverage_lca([6, 7, 8], taxa_db.parents) == 6
    assert metagenomics.coverage_lca([10, 11, 12], taxa_db.parents, 50) == 7
    assert metagenomics.coverage_lca([9], taxa_db.parents) is None
