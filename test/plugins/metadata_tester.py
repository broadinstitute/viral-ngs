#!/usr/bin/env python

"""pytest plugin for testing metadata recording, and detecting unexpected changes in commands run by each test.

Currently implemented not as a pytest plugin but as a collection of autouse fixtures.  The fixtures turn on metadata
recording for each test, and check that the metadata of commands run during the test matches the expected values.

This serves two purposes.  First, this checks that the metadata recording machinery does not interfere with normal
execution, and that the recorded metadata is reasonably correct.  Second, this serves as an additional check
against unexpected effects of code changes, that may have been missed by the tests' assertions.
This check cannot check that the output of commands run during each test is correct, but _can_ check that the
inputs, arguments and outputs of the commands have not changed.

To gather the expected command metadata for each test, we have run the test suite under various conditions while
recording command metadata for each test, then gathered the expected metadata for each test.  
See gather_expected_command_metadata_for_each_test() below.  The expected command metadata for each test is in 
test/input/TestMetadataRecording/cmd .  This data will need to be regenerated when a code change causes
the expected command metadata for some test(s) to change.
"""

import os
import os.path
import functools
import operator
import uuid
import collections
import logging

import util.misc
import util.file
import util._metadata.metadata_db as metadata_db
import util._metadata.md_utils as md_utils

import pytest

_log = logging.getLogger(__name__)

#################################################
# Gathering the expected metadata for each test
#################################################
                        
@pytest.fixture(autouse='true')
def set_nodeid_metadata(request):
    """Add pytest nodeid and a unique ID of this test invocation to the metadata"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_nodeid', request.node.nodeid), \
         util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_testid', uuid.uuid4()):
        yield

def maybe_unlist(x): 
    """If `x` is a one-element list, return just that element itself, else return the full original list"""
    return x[0] if len(x)==1 else x

def nodups(vals):
    """Return a canonicalized version of a list of values, by sorting the values, removing duplicates,
    and simplifying one-element lists to the sole element itself"""
    return maybe_unlist(list(util.misc.unique_justseen(sorted(vals, key=str))))

def unified_metadata_from_test(step_records):
    """Compute a unified representation of metadata from command(s) executed during one test.
    Records are flattened, combined into one set, and values from records with same key are merged.
    (There can be records with the same key if the test runs several commands with the same command name.)
    """
    key2vals = collections.defaultdict(list)
    for rec in step_records:
        pfx = (rec['step']['cmd_name'],)
        for k, v in util.misc.flatten_dict(rec, as_dict=(tuple, list)).items():
            key2vals[pfx+k].append(v)

    return {str(k): nodups(vals) for k, vals in key2vals.items()}

#
# Canonicalization of expected command metadata:
#
# When gathering, for regtesting purposes, the expected command metadata of commands run during a test,
# we want to factor out metadata that changes from run to run.  We do this by taking metadata for the same test
# run several times under different environments, and looking for metadata that never varies.
# Additionally, to check correctness/completeness of metadata recording, for metadata that does vary we want
# to still check that it was recorded and that the recorded values are of the correct type; we do that by
# canonicalizing the values to a default value of the value type.  For example, if a given metadata key is
# always recorded during a given test, and the value is always a string but the string is different for each run,
# we canonicalize the string value to ''.
#
# Each canonicalizer routine below canonicalizes a list of values to a canonical value; we try progressively more
# aggressive canonicalizers until we find one that maps the value of a given metadata key from all known runs of a
# given test, to the same consistent value.

def canon_identity(vals): 
    """No-op canonicalizer; applies when the set of values of a given metadata key during a given test is expected
    to always stay the same."""
    return vals

def canon_to_types(vals):
    """Map each value to a default value of the type: all strings to empty string, all ints to zero, etc."""
    return nodups([None if v is None else type(v)() for v in util.misc.make_seq(vals)])

def canon_to_None(vals):
    """If all else fails, just collapse all possible value sets to None.  This might happen if the metadata key has
    values of different types across different runs.  We'll still be checking that the given metadata
    key is expected to be present, but not checking anything about the value."""
    return None

canonicalizers = (canon_identity, canon_to_types, canon_to_None)

def regtest_fname(nodeid):
    """Return the name of the file in which the expected metadata for a given test is stored."""
    return os.path.join(util.file.get_test_input_path(), 'TestMetadataRecording', 'cmd',
                        util.file.string_to_file_name(nodeid) + '.json.gz')

def gather_expected_command_metadata_for_each_test():
    """For each test, gather the expected command metadata for commands run during that test, for regtesting purposes.
    
    This function is called when this module is run as a top-level script.  Before running it, run the test suite
    under several conditions while recording metadata, so that for each test we have the metadata record of commands
    run by that test, and can detect which parts of metadata vary and which can be expected to stay the same for 
    future tests.  We save the expected metadata for each test to test/input/TestMetadataRecording/cmd .
    We err on the side of avoiding false alarms, rather than on the side of detecting all minute metadata changes.
    """

    # First, gather all comand metadata records that we have, and group them by the pytest node id, and within that,
    # group commands executed within a single invocation of a that node id test.
    
    recs = [rec for rec in metadata_db.load_all_records() if md_utils.dict_has_keys(rec['step']['metadata_from_cmd_line'],
                                                                                    'pytest_nodeid pytest_testid')]
    nodeid2testid2recs = collections.defaultdict(functools.partial(collections.defaultdict, list))
    for rec in recs:
        mdata = rec['step']['metadata_from_cmd_line']
        nodeid2testid2recs[mdata['pytest_nodeid']][mdata['pytest_testid']].append(rec)

    # Now identify for each node ID (1) what metadata keys are always present and (2) of those, which always have
    # the same set of values for a test, and which have a set of values that can vary from run to run.
    # For a given metadadata key, there can be multiple values for a given test if that test runs multiple
    # instances of a command with a given name.

    for nodeid, testid2recs in nodeid2testid2recs.items():
        unified_recs = list(map(unified_metadata_from_test, testid2recs.values()))
        _log.info('for nodeid {} have {} recs'.format(nodeid, len(unified_recs)))
        if len(unified_recs) < 2: 
            _log.info('skipping nodeid {} because too few recs'.format(nodeid))
            continue

        # Only check metadata keys seen in all known runs of the given nodeid.
        # If later code changes add new keys we won't report an error, but if later code changes drop 
        # of the previously expected keys we will.
        common_keys = set.intersection(*[set(k2v.keys()) for k2v in unified_recs])
        canonicalizer2key2vals = collections.defaultdict(functools.partial(collections.defaultdict, list))
        for k in common_keys:
            # across all known runs of this nodeid, what values of this metadata key have we seen?
            vals = [ k2v[k] for k2v in unified_recs ]

            # find a canonicalizer under which, in all executions of this nodeid, the set of values
            # for metadata key k is the same.  Try increasingly aggressive canonicalizers.
            for canonicalizer in canonicalizers:
                canon_vals = list(map(canonicalizer, vals))
                if all(v == canon_vals[0] for v in canon_vals):
                    canonicalizer2key2vals[canonicalizer.__name__][k] = canon_vals[0]
                    break

        if canonicalizer2key2vals:
            fname = regtest_fname(nodeid)
            with open(fname, 'wb') as out:
                out.write(util.file.to_json_gz(canonicalizer2key2vals, filename=fname))

####################################################################
# Checking of command metadata from each test against expectation
####################################################################

#@pytest.fixture(autouse='true')
@pytest.fixture()
def check_metadata_of_cmds_in_test(request, tmpdir_factory):
    """Checks that the metadata of commands executed during each test has not changed.
    
    For each test, we have created and stored under test/input/TestMetadataRecording/cmd a record of metadata
    expected to be produced by commands run during that test.  (Of course, only for tests that run one or more commands.)
    """
    metadata_db_path = tmpdir_factory.mktemp('metadata_db')

    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path, sep=' '):
        yield metadata_db_path
        with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):  # load just the metadata from this test
            recs = unified_metadata_from_test(metadata_db.load_all_records())
            regtest_file = regtest_fname(request.node.nodeid)
            if os.path.isfile(regtest_file):
                with open(regtest_file, 'rb') as rf:
                    regtest_data = util.file.from_json_gz(rf.read())
                for canonicalizer, key2vals in regtest_data.items():
                    canonicalizer_fn = globals()[canonicalizer]
                    for k, v in key2vals.items():
                        msg=('Something changed about the set of commands run during test {}: the commands themselves, '
                             'their inputs or outputs, or their parameters.  Metadata for key {} is different from '
                             'expected.  If this change is intended, please update the '
                             'expected metadata in {}, '
                             'by running test/plugins/metadata_tester.py, or just delete the expected metadata file '
                             'to stop testing metadata consistency for this test.'.format(request.node.nodeid, k, 
                                                                                          regtest_file))
                        assert k in recs, msg
                        assert canonicalizer_fn(recs[k]) == v, msg

###################################################################                

if __name__ == '__main__':
    gather_expected_command_metadata_for_each_test()

