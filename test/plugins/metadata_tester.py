"""pytest plugin for testing metadata recording"""

import os
import os.path
import functools
import operator
import uuid
import collections

import util.misc
import util.file
import util._metadata.metadata_db as metadata_db
import util._metadata.md_utils as md_utils

import pytest

def canonicalize_step_record(step_record, canonicalize_keys=()):
    """Return a canonicalized flat dict of key-value pairs representing this record, for regression testing purposes.
    Canonicalization uniformizes or drops fields that may change between runs.
    """
    def canonicalize_value(v):
        if v is None: return v
        return type(v)()
    pfx = (step_record['step']['cmd_name'],)
    return {pfx+k: canonicalize_value(v) \
            if md_utils.tuple_key_matches(k, canonicalize_keys) \
            else v for k, v in util.misc.flatten_dict(step_record, as_dict=(tuple,list)).items() \
            if k[:3] != ('step', 'run_info', 'argv')}

def canonicalize_step_records(step_records, canonicalize_keys=()):
    """Canonicalize a group of step records"""
    return sorted(map(str, functools.reduce(operator.concat, 
                                            [list(canonicalize_step_record(r, canonicalize_keys).items())
                                             for r in step_records], [])))

@pytest.fixture(scope='session', autouse='true')
def tmp_metadata_db(tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = os.environ.get('VIRAL_NGS_TEST_METADATA_PATH', tmpdir_factory.mktemp('metadata_db'))
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path, sep=' '):
        yield metadata_db_path        


#@pytest.fixture(autouse='true')
def per_test_metadata_db(request, tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir_factory.mktemp('metadata_db')

    # The following tests, for varying reasons, result in some captured metadata varying between test runs.
    # For now, just don't attempt to check that the metadata does not change; eventually some/most of these exceptions
    # can be removed.
    nondet_tests = (
        'test/unit/test_assembly.py::TestDeambigAndTrimFasta::test_deambig_fasta',
    )

    canonicalize_keys = (
        ('step', 'run_env run_info run_id step_id version_info'),
        ('step', 'args', '', 'val'),
        ('step', 'args', '', ' '.join(map(str, range(10))), 'val'),
        ('step', 'args', '', 'files', '',
         'abspath ctime device fname inode mtime owner realpath'),
        ('step', 'args', '', ' '.join(map(str, range(10))), 'files', '',
         'abspath ctime device fname inode mtime owner realpath'),
        ('step', 'metadata_from_cmd_return', 'runtime'),
        ('step', 'enclosing_steps'),
        ('step', 'args', 'tmp_dir'),
        ('step', 'args', 'tmp_dirKeep'),
        ('step', 'args', 'novo_params'),
        ('step', 'args', 'refDbs'),
        ('step', 'args', 'outputDirectory', 'files', '2'),
        ('step', 'args', 'inVcf', 'files', '0'),
        ('step', 'args', 'blastDbs', '0', 'files', '2'),
    )

    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path, sep=' '):
        print('metadata_db_path:', metadata_db_path)
        yield metadata_db_path
        with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
            recs = unified_metadata_from_test(metadata_db.load_all_records())
            regtest_file = regtest_fname(request.node.nodeid)
            if os.path.isfile(regtest_file):
                with open(regtest_file, 'rb') as rf:
                    regtest_data = util.file.from_json_gz(rf.read())
                for canonicalizer, key2vals in regtest_data.items():
                    canonicalizer_fn = globals()[canonicalizer]
                    for k, v in key2vals.items():
                        assert k in recs
                        assert canonicalizer_fn(recs[k]) == v

@pytest.fixture(scope='session', autouse='true')
def no_detailed_env():
    """Disable time-consuming gathering of detailed env"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_DETAILED_ENV', None):
        yield

@pytest.fixture(autouse='true')
def set_nodeid_metadata(request):
    """Add pytest nodeid and a unique ID of this test invocation to the metadata"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_nodeid', request.node.nodeid), \
         util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_testid', uuid.uuid4()):
        yield

def maybe_unlist(x): return x[0] if len(x)==1 else x

def unified_metadata_from_test(step_records):
    """Compute a unified representation of metadata from command(s) executed during one test.
    Records are flattened, combined into one set, and values from records with same key are merged."""
    key2vals = collections.defaultdict(list)
    for rec in step_records:
        pfx = (rec['step']['cmd_name'],)
        for k, v in util.misc.flatten_dict(rec, as_dict=(tuple, list)).items():
            key2vals[pfx+k].append(v)

    return {str(k): maybe_unlist(list(util.misc.unique_justseen(sorted(vals, key=str)))) for k, vals in key2vals.items()}

def canon_identity(vals): return vals

def canon_to_types(vals):
    return maybe_unlist([None if v is None else type(v)() for v in util.misc.make_seq(vals)])

def canon_to_None(vals):
    return None

canonicalizers = (canon_identity, canon_to_types, canon_to_None)

def regtest_fname(nodeid):
    return os.path.join(util.file.get_test_input_path(), 'TestMetadataRecording', 'cmd',
                        util.file.string_to_file_name(nodeid) + '.json.gz')

def gather_metadata_regtests():
    """Gather data for metadata regtesting."""
    
    recs = [rec for rec in metadata_db.load_all_records() if md_utils.dict_has_keys(rec['step']['metadata_from_cmd_line'],
                                                                                    'pytest_nodeid pytest_testid')]
    nodeid2testid2recs = collections.defaultdict(functools.partial(collections.defaultdict, list))
    for rec in recs:
        mdata = rec['step']['metadata_from_cmd_line']
        nodeid2testid2recs[mdata['pytest_nodeid']][mdata['pytest_testid']].append(rec)

    for nodeid, testid2recs in nodeid2testid2recs.items():
        unified_recs = list(map(unified_metadata_from_test, testid2recs.values()))
        print('for nodeid {} have {} recs'.format(nodeid, len(unified_recs)))
        if len(unified_recs) < 2: continue
        common_keys = set.intersection(*[set(k2v.keys()) for k2v in unified_recs])
        canonicalizer2key2vals = collections.defaultdict(functools.partial(collections.defaultdict, list))
        for k in common_keys:
            vals = [ k2v[k] for k2v in unified_recs ]
            for canonicalizer in canonicalizers:
                canon_vals = list(map(canonicalizer, vals))
                if all(v == canon_vals[0] for v in canon_vals):
                    canonicalizer2key2vals[canonicalizer.__name__][k] = canon_vals[0]
                    break

        if canonicalizer2key2vals:
            fname = regtest_fname(nodeid)
            with open(fname, 'wb') as out:
                out.write(util.file.to_json_gz(canonicalizer2key2vals, filename=fname))

if __name__ == '__main__':
    gather_metadata_regtests()
