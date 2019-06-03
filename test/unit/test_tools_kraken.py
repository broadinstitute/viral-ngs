# Unit tests for kraken
import os.path

import pytest

import util.file
import util.misc
import tools.kraken
from test import _CPUS


@pytest.fixture
def krakenuniq():
    return tools.kraken.KrakenUniq()


@pytest.fixture
def in_bam():
    return os.path.join(util.file.get_test_input_path(), 'almost-empty.bam')


@pytest.fixture
def db(tmpdir_factory):
    return str(tmpdir_factory.mktemp('db'))


@pytest.fixture(autouse=True)
def mocks(mocker):
    mock_run = mocker.patch('util.misc.run', autospec=True)
    mock_check_call = mocker.patch('subprocess.check_call', autospec=True)
    return {
        'run': mock_run,
        'check_call': mock_check_call,
    }


def test_krakenuniq_classify(mocks, krakenuniq, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')
    krakenuniq.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'krakenuniq' == os.path.basename(args[0])
    assert util.misc.list_contains(['--db', db], args)
    assert util.misc.list_contains(['--output', out_reads], args)
    assert util.misc.list_contains(['--threads', str(_CPUS)], args)

def test_classify_num_threads(mocks, krakenuniq, db, in_bam):
    out_reads = util.file.mkstempfname('.reads.txt')

    krakenuniq.classify(in_bam, db, out_reads)
    args = mocks['check_call'].call_args[0][0]
    assert 'krakenuniq' == os.path.basename(args[0])
    assert '--threads' in args
    actual = args[args.index('--threads')+1]
    assert actual == str(_CPUS)

    for requested in (1,2,3,8,11,20):
        expected = min(_CPUS, requested)
        krakenuniq.classify(in_bam, db, out_reads, num_threads=requested)
        args = mocks['check_call'].call_args[0][0]
        assert 'krakenuniq' == os.path.basename(args[0])
        assert '--threads' in args
        actual = args[args.index('--threads')+1]
        assert actual == str(expected), "failure for requested %s, expected %s, actual %s" % (requested, expected, actual)
