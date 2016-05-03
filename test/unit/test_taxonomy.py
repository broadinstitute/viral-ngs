from collections import Counter
import io
import os.path
import textwrap
import pytest
import util.file
import taxonomy


@pytest.fixture
def taxa_db_simple():
    db = taxonomy.TaxonomyDb()
    db.gis = {1:2, 2:3, 3:4, 4:5}
    db.parents = {1: 1, 2: 1, 3: 2, 4: 3, 5: 4}
    return db


@pytest.fixture
def taxa_db(parents, names, ranks):
    db = taxonomy.TaxonomyDb()
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
    test_path = os.path.join(util.file.get_test_input_path(),
                             'TestTaxonomy')
    return open(os.path.join(test_path, 'simple.m8'))


def test_tree_level_lookup(parents):
    level_cache = {1: 1}
    assert taxonomy.tree_level_lookup(parents, 1, level_cache) == 1
    assert taxonomy.tree_level_lookup(parents, 3, level_cache) == 2
    assert taxonomy.tree_level_lookup(parents, 12, level_cache) == 6
    level_cache = {1: 1}
    assert taxonomy.tree_level_lookup(parents, 12, level_cache) == 6
    assert taxonomy.tree_level_lookup(parents, 8, level_cache) == 5


def test_push_up_tree_hits(parents):
    hits = Counter({1: 3, 3: 5, 6: 3, 7: 3, 13: 5})
    with pytest.raises(AssertionError):
        taxonomy.push_up_tree_hits(parents, hits)

    expected = hits.copy()
    assert taxonomy.push_up_tree_hits(parents, hits.copy(), min_support=1) == expected

    expected = Counter({3: 5, 6: 6, 13: 5})
    assert taxonomy.push_up_tree_hits(parents, hits.copy(), min_support=5) == expected

    assert (taxonomy.push_up_tree_hits(parents, hits.copy(), min_support=10) ==
            Counter({6: 11}))
    assert (taxonomy.push_up_tree_hits(parents, hits.copy(), min_support=18) ==
            Counter({1: 19}))
    assert (taxonomy.push_up_tree_hits(parents, hits.copy(), min_support_percent=100) ==
            Counter({1: 19}))


def test_parents_to_children(parents):
    children = taxonomy.parents_to_children(parents)
    assert children[1] == [3]


def test_rank_code():
    assert taxonomy.rank_code('species') == 'S'
    assert taxonomy.rank_code('genus') == 'G'
    assert taxonomy.rank_code('superkingdom') == 'D'
    assert taxonomy.rank_code('nonexist') == '-'


def test_blast_records(simple_m8):
    test_path = os.path.join(util.file.get_test_input_path(),
                             'TestTaxonomy')
    with simple_m8 as f:
        records = list(taxonomy.blast_records(f))
    assert len(records) == 110
    assert records[0].bit_score == 63.5
    assert records[-1].bit_score == 67.4


def test_blast_lca(taxa_db_simple, simple_m8):
    test_path = os.path.join(util.file.get_test_input_path(),
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
    out = io.StringIO()
    with simple_m8 as f:
        taxonomy.blast_lca(taxa_db_simple, f, out, paired=True)
        out.seek(0)
        assert out.read() == expected


def test_paired_query_id():
    tup = ['query', 'gi|10|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5]

    blast1 = taxonomy.BlastRecord(*tup)
    assert taxonomy.paired_query_id(blast1) == blast1

    new_tup = tup.copy()
    new_tup[0] = 'query/1'
    new_blast1 = taxonomy.BlastRecord(*new_tup)
    assert taxonomy.paired_query_id(new_blast1) == blast1

    new_tup = tup.copy()
    new_tup[0] = 'query/2'
    new_blast1 = taxonomy.BlastRecord(*new_tup)
    assert taxonomy.paired_query_id(new_blast1) == blast1

    new_tup = tup.copy()
    new_tup[0] = 'query/3'
    new_blast1 = taxonomy.BlastRecord(*new_tup)
    assert taxonomy.paired_query_id(new_blast1) == new_blast1


def test_translate_gi_to_tax_id(taxa_db_simple):
    tup = ['query', 'gi|4|else', 90., 80, 60, 2, 30, 80,
           1100, 1150, 1e-7, 64.5]
    blast1 = taxonomy.BlastRecord(*tup)

    tup[1] = 5
    expected = taxonomy.BlastRecord(*tup)
    assert taxonomy.translate_gi_to_tax_id(taxa_db_simple, blast1) == expected


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
    report = taxonomy.kraken_dfs_report(taxa_db, hits)
    text_report = '\n'.join(list(report)) + '\n'
    assert text_report == expected
