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

from unittest import mock
from unittest.mock import patch

from viral_ngs.core import picard
from viral_ngs import metagenomics
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc
from tests import TestCaseWithTmp, assert_equal_bam_reads, _CPUS

if six.PY2:
    from StringIO import StringIO
else:
    from io import StringIO


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in metagenomics.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestKronaCalls(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = patch('viral_ngs.metagenomics.krona.Krona', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_krona = patcher.start()

        self.inTsv = util_file.mkstempfname('.tsv')
        self.db = tempfile.mkdtemp('db')

    def test_krona_import_taxonomy(self):
        out_html = util_file.mkstempfname('.html')
        metagenomics.main_krona([self.inTsv], self.db, out_html, queryColumn=3, taxidColumn=5, scoreColumn=7,
                           noHits=True, noRank=True, inputType='tsv')
        self.mock_krona().import_taxonomy.assert_called_once_with(
            self.db, [self.inTsv + ',' + os.path.basename(self.inTsv)], out_html, query_column=3, taxid_column=5, score_column=7,
            no_hits=True, no_rank=True, magnitude_column=None)


def test_centrifuger_parser_invokes_tool():
    with patch('viral_ngs.metagenomics.centrifuger.Centrifuger', autospec=True) as mock_centrifuger:
        args = [
            'db',
            'input.bam',
            'out.tsv',
            '--k', '3',
            '--unclassified_prefix', 'unclassified',
            '--classified_prefix', 'classified',
            '--min_hitlen', '17',
            '--hitk_factor', '5',
            '--merge_readpair',
            '--threads', '4',
        ]
        args = metagenomics.parser_centrifuger(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_centrifuger.return_value.classify.assert_called_once_with(
            'input.bam',
            'db',
            'out.tsv',
            k=3,
            unclassified_prefix='unclassified',
            classified_prefix='classified',
            min_hitlen=17,
            hitk_factor=5,
            merge_readpair=True,
            num_threads=4,
        )


def test_centrifuger_parser_defaults_to_k_one():
    with patch('viral_ngs.metagenomics.centrifuger.Centrifuger', autospec=True) as mock_centrifuger:
        args = metagenomics.parser_centrifuger(argparse.ArgumentParser()).parse_args([
            'db',
            'input.bam',
            'out.tsv',
        ])
        args.func_main(args)

        mock_centrifuger.return_value.classify.assert_called_once_with(
            'input.bam',
            'db',
            'out.tsv',
            k=1,
            unclassified_prefix=None,
            classified_prefix=None,
            min_hitlen=None,
            hitk_factor=None,
            merge_readpair=False,
            num_threads=_CPUS,
        )


def test_centrifuger_build_parser_invokes_tool():
    with patch('viral_ngs.metagenomics.centrifuger.Centrifuger', autospec=True) as mock_centrifuger:
        args = [
            'db_prefix',
            'nodes.dmp',
            'names.dmp',
            '--ref_fastas', 'ref1.fna', 'ref2.fna',
            '--conversion_table', 'seqid_to_taxid.map',
            '--build_mem', '256M',
            '--threads', '4',
        ]
        args = metagenomics.parser_centrifuger_build(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_centrifuger.return_value.build.assert_called_once_with(
            'db_prefix',
            'nodes.dmp',
            'names.dmp',
            ref_fastas=['ref1.fna', 'ref2.fna'],
            ref_list=None,
            conversion_table='seqid_to_taxid.map',
            build_mem='256M',
            num_threads=4,
        )


def test_centrifuger_quant_parser_invokes_tool():
    with patch('viral_ngs.metagenomics.centrifuger.Centrifuger', autospec=True) as mock_centrifuger:
        args = [
            'db',
            'classification.tsv',
            'quant.tsv',
            '--min_score', '10',
            '--min_length', '50',
            '--output_format', '1',
        ]
        args = metagenomics.parser_centrifuger_quant(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_centrifuger.return_value.quant.assert_called_once_with(
            'db',
            'classification.tsv',
            'quant.tsv',
            min_score=10,
            min_length=50,
            output_format=1,
        )


def test_centrifuger_classification_to_kraken2_parser_invokes_tool():
    with patch(
            'viral_ngs.metagenomics.centrifuger.Centrifuger.classification_to_kraken2',
            autospec=True,
    ) as mock_classification_to_kraken2:
        args = [
            'classification.tsv',
            'kraken2.tsv',
        ]
        args = metagenomics.parser_centrifuger_classification_to_kraken2(
            argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_classification_to_kraken2.assert_called_once_with(
            'classification.tsv',
            'kraken2.tsv',
        )


def test_centrifuger_kreport_parser_invokes_tool():
    with patch('viral_ngs.metagenomics.centrifuger.Centrifuger', autospec=True) as mock_centrifuger:
        args = [
            'db',
            'classification.tsv',
            'kreport.tsv',
            '--no_lca',
            '--show_zeros',
            '--is_count_table',
            '--min_score', '10',
            '--min_length', '50',
            '--report_score_data',
        ]
        args = metagenomics.parser_centrifuger_kreport(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_centrifuger.return_value.kreport.assert_called_once_with(
            'db',
            'classification.tsv',
            'kreport.tsv',
            no_lca=True,
            show_zeros=True,
            is_count_table=True,
            min_score=10,
            min_length=50,
            report_score_data=True,
        )


class TestVirNucProCalls(TestCaseWithTmp):
    @patch('viral_ngs.metagenomics.virnucpro.classify_contigs', autospec=True)
    def test_virnucpro_contigs(self, mock_classify_contigs):
        args = metagenomics.parser_virnucpro_contigs(
            argparse.ArgumentParser()
        ).parse_args([
            'highestscore.tsv',
            'contigs.tsv',
            '--min-viral-prop', '0.2',
            '--min-nonviral-prop', '0.3',
            '--min-chunks', '7',
            '--min-confident-score', '0.75',
            '--max-opposing-score', '0.25',
            '--min-ambiguous-score', '0.65',
            '--min-weighted-delta', '0.35',
            '--high-confidence-delta', '0.55',
            '--id-col', 'Modified_ID',
            '--id-pattern', r'(NODE_\d+)',
        ])
        args.func_main(args)

        mock_classify_contigs.assert_called_once_with(
            'highestscore.tsv',
            'contigs.tsv',
            min_viral_prop=0.2,
            min_nonviral_prop=0.3,
            min_chunks=7,
            min_confident_score=0.75,
            max_opposing_score=0.25,
            min_ambiguous_score=0.65,
            min_weighted_delta=0.35,
            high_confidence_delta=0.55,
            id_col='Modified_ID',
            id_pattern=r'(NODE_\d+)',
        )

    @patch('viral_ngs.metagenomics.virnucpro.classify_reads_by_contig', autospec=True)
    def test_virnucpro_label_reads_by_contig(self, mock_classify_reads_by_contig):
        args = metagenomics.parser_virnucpro_label_reads_by_contig(
            argparse.ArgumentParser()
        ).parse_args([
            'reads.bam',
            'contigs.tsv',
            'reads_classified.tsv.zst',
            '--min-mapq', '10',
            '--min-identity', '95.0',
            '--min-query-cov', '85.0',
            '--duckdb-memory-limit', '4GB',
            '--work-dir', '/tmp/virnucpro',
        ])
        args.func_main(args)

        mock_classify_reads_by_contig.assert_called_once_with(
            'reads.bam',
            'contigs.tsv',
            'reads_classified.tsv.zst',
            min_mapq=10,
            min_identity=95.0,
            min_query_cov=85.0,
            duckdb_memory_limit='4GB',
            work_dir='/tmp/virnucpro',
        )

    def test_virnucpro_label_reads_by_contig_help_uses_percent_units(self):
        help_text = metagenomics.parser_virnucpro_label_reads_by_contig(
            argparse.ArgumentParser()).format_help()

        assert 'Fractional-scale' not in help_text
        assert 'percent units, e.g. 90 not 0.9' in help_text
        assert 'percent units, e.g. 80 not 0.8' in help_text


def test_kallisto_top_taxa_parser_invokes_tool():
    with patch('viral_ngs.metagenomics.kallisto.Kallisto', autospec=True) as mock_kallisto:
        args = [
            'counts.tsv',
            'top_taxa.tsv',
            '--id-to-tax-map', 'taxonomy.csv',
            '--target-taxon', 'Viruses',
        ]
        args = metagenomics.parser_kallisto_top_taxa(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        mock_kallisto.return_value.write_top_taxa_report_from_counts_tsv.assert_called_once_with(
            'counts.tsv',
            'top_taxa.tsv',
            id_to_tax_map='taxonomy.csv',
            target_taxon='Viruses',
        )


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
    test_path = join(util_file.get_test_input_path(),
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
    test_path = join(util_file.get_test_input_path(),
                             'TestTaxonomy')
    with simple_m8 as f:
        records = list(metagenomics.blast_records(f))
    assert len(records) == 110
    assert records[0].bit_score == 63.5
    assert records[-1].bit_score == 67.4


def test_blast_lca(taxa_db_simple, simple_m8):
    test_path = join(util_file.get_test_input_path(),
                             'TestTaxonomy')
    expected = textwrap.dedent("""\
    C\tM04004:13:000000000-AGV3H:1:1101:12068:2105\t2
    C\tM04004:13:000000000-AGV3H:1:1101:13451:2146\t2
    C\tM04004:13:000000000-AGV3H:1:1101:13509:2113\t2
    C\tM04004:13:000000000-AGV3H:1:1101:14644:2160\t2
    C\tM04004:13:000000000-AGV3H:1:1101:18179:2130\t2
    C\tM04004:13:000000000-AGV3H:1:1111:10629:2610\t2
    C\tM04004:13:000000000-AGV3H:1:1111:10629:26101\t2
    """)
    out = StringIO()
    with simple_m8 as f:
        taxa_db_simple.blast_lca(f, out, paired=True)
        out.seek(0)
        assert out.read() == expected


def test_paired_query_id():
    tup = ['query', 'gi|10|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5, []]

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
           1100, 1150, 1e-7, 64.5, []]
    blast1 = metagenomics.BlastRecord(*tup)

    tup[1] = 5
    expected = metagenomics.BlastRecord(*tup)
    assert taxa_db_simple.translate_gi_to_tax_id(blast1) == expected


def test_ancestor_lookup(taxa_db_simple):
    assert taxa_db_simple.get_ordered_ancestors(4) == [3, 2, 1]


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
    report = taxa_db.kraken_dfs_report(hits)
    text_report = '\n'.join(list(report)) + '\n'
    assert text_report == expected


def test_coverage_lca(taxa_db):
    assert metagenomics.coverage_lca([10, 11, 12], taxa_db.parents) == 6
    assert metagenomics.coverage_lca([1, 3], taxa_db.parents) == 1
    assert metagenomics.coverage_lca([6, 7, 8], taxa_db.parents) == 6
    assert metagenomics.coverage_lca([10, 11, 12], taxa_db.parents, 50) == 7
    assert metagenomics.coverage_lca([9], taxa_db.parents) is None

@pytest.mark.skip(reason="KrakenUniq disabled from future versions for now, pending conda rebuild of @yesimon's custom fork")
def test_krakenuniq():
    with patch('viral_ngs.metagenomics.kraken.KrakenUniq.pipeline') as p:
        args = [
            'db',
            'input.bam',
            '--outReports', 'output.report',
            '--outReads', 'output.reads',
        ]
        args = metagenomics.parser_krakenuniq(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)
        p.assert_called_with('db', ['input.bam'], num_threads=mock.ANY, filter_threshold=mock.ANY, out_reports=['output.report'], out_reads=['output.reads'])


class TestBamFilter(TestCaseWithTmp):
    def test_bam_filter_simple(self):
        temp_dir = tempfile.gettempdir()
        input_dir = util_file.get_test_input_path(self)
        taxonomy_dir = os.path.join(util_file.get_test_input_path(),"TestMetagenomicsSimple","db","taxonomy")

        filtered_bam = util_file.mkstempfname('.bam')
        args = [
            os.path.join(input_dir,"input.bam"),
            os.path.join(input_dir,"input.kraken-reads.tsv.gz"),
            filtered_bam,
            os.path.join(taxonomy_dir,"nodes.dmp"),
            os.path.join(taxonomy_dir,"names.dmp"),
            "--taxNames",
            "Ebolavirus"
        ]
        args = metagenomics.parser_filter_bam_to_taxa(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        expected_bam = os.path.join(input_dir,"expected.bam")
        assert_equal_bam_reads(self, filtered_bam, expected_bam)

    def test_bam_filter_by_tax_id(self):
        temp_dir = tempfile.gettempdir()
        input_dir = util_file.get_test_input_path(self)
        taxonomy_dir = os.path.join(util_file.get_test_input_path(),"TestMetagenomicsSimple","db","taxonomy")

        filtered_bam = util_file.mkstempfname('.bam')
        args = [
            os.path.join(input_dir,"input.bam"),
            os.path.join(input_dir,"input.kraken-reads.tsv.gz"),
            filtered_bam,
            os.path.join(taxonomy_dir,"nodes.dmp"),
            os.path.join(taxonomy_dir,"names.dmp"),
            "--taxIDs",
            "186538"
        ]
        args = metagenomics.parser_filter_bam_to_taxa(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        expected_bam = os.path.join(input_dir,"expected.bam")
        assert_equal_bam_reads(self, filtered_bam, expected_bam)
