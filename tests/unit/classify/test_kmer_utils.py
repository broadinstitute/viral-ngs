"""Unit tests for kmer_utils.py"""

__author__ = "ilya@broadinstitute.org"

import os
import sys
import collections
import argparse
import logging
import itertools
import unittest

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from viral_ngs import kmer_utils
from viral_ngs import read_utils
from viral_ngs.core import cmd as util_cmd
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc
from viral_ngs.classify import kmc
from viral_ngs.core import samtools


_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in kmer_utils.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

#################################
# Some general utils used below #
#################################

def _seq_as_str(s):  # pylint: disable=invalid-name
    """Return a sequence as a str, regardless of whether it was a str, a Seq or a SeqRecord"""
    if isinstance(s, Seq):
        return str(s)
    if isinstance(s, SeqRecord):
        return str(s.seq)
    return s

def _yield_seq_recs(seq_file):
    """Yield sequence records from the file, regardless of file format."""
    with util_file.tmp_dir(suffix='_seqs_as_strs') as t_dir:
        if seq_file.endswith('.bam'):
            t_fa = os.path.join(t_dir, 'bam2fa.fasta')
            samtools.SamtoolsTool().bam2fa(seq_file, t_fa, append_mate_num=True)
            seq_file = t_fa
        with util_file.open_or_gzopen(seq_file, 'rt') as seq_f:
            for rec in Bio.SeqIO.parse(seq_f, util_file.uncompressed_file_type(seq_file)[1:]):
                yield rec

def _list_seq_recs(seq_file):
    """Return a list of sequence records from the file, regardless of file format."""
    return list(_yield_seq_recs(seq_file))


def _yield_seqs_as_strs(seqs):
    """Yield sequence(s) from `seqs` as strs.  seqs can be a str/SeqRecord/Seq, a filename of a sequence file,
    or an iterable of these.  If a filename of a sequence file, all sequences from that file are yielded."""
    for seq in util_misc.make_seq(seqs, (str, SeqRecord, Seq)):
        seq = _seq_as_str(seq)
        if not any(seq.endswith(ext) for ext in '.fasta .fasta.gz .fastq .fastq.gz .bam'):
            yield seq
        else:
            for rec in _yield_seq_recs(seq):
                yield str(rec.seq)

def _list_seqs_as_strs(seqs):
    """Return a list of sequence(s) from `seqs` as strs.  seqs can be a str/SeqRecord/Seq, a filename of a sequence file,
    or an iterable of these.  If a filename of a sequence file, all sequences from that file are yielded."""
    return list(_yield_seqs_as_strs(seqs))

def _getargs(args, valid_args):
    """Extract valid args from an argparse.Namespace or a dict.  Returns a dict containing keys from `args`
    that are in `valid_args`; `valid_args` is a space-separated list of valid args."""
    return util_misc.subdict(vars(args) if isinstance(args, argparse.Namespace) else args, valid_args.split())

def _strip_mate_num(rec_id):
    """Name of a read, with any /1 or /2 at the end removed"""
    return rec_id[:-2] if rec_id.endswith('/1') or rec_id.endswith('/2') else rec_id

#################################################################################################################


class KmcPy(object):
    """Reimplementation of some kmc functions in simple Python code.

    To help generate the expected correct output, we reimplement in simple Python code some
    of KMC's functionality.   This is also to make up for KMC's lack of a public test suite
    ( https://github.com/refresh-bio/KMC/issues/55 ).
    """

    def _revcomp(self, kmer):
        """Return the reverse complement of a kmer, given as a string"""
        assert isinstance(kmer, str)
        return str(SeqRecord(Seq(kmer), annotations={"molecule_type": "DNA"}).seq.reverse_complement())

    def _canonicalize(self, kmer):
        """Return the canonical version of a kmer"""
        return min(kmer, self._revcomp(kmer))

    def _compute_kmers_iter(self, seq_strs, kmer_size, single_strand, **ignore):
        """Yield kmers of seq(s).  Unless `single_strand` is True, each kmer
        is canonicalized before being returned.  Note that
        deduplication is not done, so each occurrence of each kmer is
        yielded.  Kmers containing non-TCGA bases are skipped.
        """
        for seq in seq_strs:
            n_kmers = len(seq)-kmer_size+1

            # mark kmers containing invalid base(s)
            valid_kmer = [True] * n_kmers
            for i in range(len(seq)):  # pylint: disable=consider-using-enumerate
                if seq[i].upper() not in 'TCGA':
                    invalid_lo = max(0, i-kmer_size+1)
                    invalid_hi = min(n_kmers, i+1)
                    valid_kmer[invalid_lo:invalid_hi] = [False]*(invalid_hi-invalid_lo)

            for i in range(n_kmers):
                if valid_kmer[i]:
                    kmer = seq[i:i+kmer_size]
                    yield kmer if single_strand else self._canonicalize(kmer)

    def _compute_kmers(self, *args, **kw):
        """Return list of kmers of seq(s).  Unless `single_strand` is True, each kmer
        is canonicalized before being returned.  Note that
        deduplication is not done, so each occurrence of each kmer is
        yielded.  Kmers containing non-TCGA bases are skipped.
        """
        return list(self._compute_kmers_iter(*args, **kw))

    def compute_kmer_counts(self, seq_files, kmer_size, min_occs, max_occs,
                            counter_cap, single_strand, **ignore):
        """Yield kmer counts of seq(s).  Unless `single_strand` is True, each kmer is
        canonicalized before being counted.  Kmers containing non-TCGA bases are skipped.
        Kmers with fewer than `min_occs` or more than `max_occs` occurrences
        are dropped, and kmer counts capped at `counter_cap`, if these args are given.
        """
        counts = collections.Counter(self._compute_kmers(_list_seqs_as_strs(seq_files), kmer_size, single_strand))
        return self._filter_kmer_counts(counts=counts, min_occs=min_occs, max_occs=max_occs, counter_cap=counter_cap)

    def _filter_kmer_counts(self, counts, min_occs=None, max_occs=None, counter_cap=None):
        """From a dict of kmer counts, drop kmers with counts below `min_occs` or above `max_occs`, and
        cap counter values at `counter_cap`."""
        return collections.Counter({kmer : min(count, counter_cap or count) \
                                    for kmer, count in counts.items() \
                                    if (count >= (min_occs or count)) and \
                                    (count <= (max_occs or count))})


    def filter_reads(self, db_kmer_counts, in_reads, kmer_size, single_strand,
                    db_min_occs, db_max_occs, read_min_occs, read_max_occs,
                    read_min_occs_frac, read_max_occs_frac, **ignore):
        """Fiter sequences based on their kmer contents; returns the ids of the passing seqs.

        Inputs:
           in_reads: sequences to filter, as a fasta/fastq/bam

        Params:
           db_kmer_counts: a collections.Counter mapping kmers to their counts in the kmer database
           kmer_size: kmer size
           single_strand: if False, kmers are canonicalized
           db_min_occs: drop from db_kmer_counts kmers with counts below this
           db_max_occs: drop from db_kmer_counts kmers with counts above this
           read_min_occs: drop reads with fewer than this many occurrences of kmers from the database
           read_max_occs: drop reads with more than this many occurrences of kmers from the database
           read_min_occs_frac: only keep reads with at least this many occurrences of kmers from database,
             interpreted as a fraction of read length in kmers
           read_max_occs_frac: only keep reads with no more than this many occurrence of kmers from the database.
             interpreted as a fraction of read length in kmers.
        """
        db_kmers = self._filter_kmer_counts(counts=db_kmer_counts, min_occs=db_min_occs, max_occs=db_max_occs).keys()

        seqs_ids_out = set()
        rel_thresholds = (read_min_occs_frac, read_max_occs_frac) != (0., 1.)
        in_recs = _list_seq_recs(in_reads)

        seq_occs_hist = collections.Counter()
        mate_cnt = collections.Counter()

        for rec in in_recs:
            seq = str(rec.seq)
            seq_kmer_counts = self.compute_kmer_counts(seq, kmer_size=kmer_size, single_strand=single_strand,
                                                       min_occs=1, max_occs=util_misc.MAX_INT32,
                                                       counter_cap=util_misc.MAX_INT32)
            assert not single_strand
            seq_occs = sum([seq_count for kmer, seq_count in seq_kmer_counts.items() \
                            if kmer in db_kmers])
            seq_occs_hist[seq_occs] += 1

            if rel_thresholds:
                n_seq_kmers = len(seq)-kmer_size+1
                read_min_occs_seq, read_max_occs_seq = (int(read_min_occs_frac * n_seq_kmers),
                                                        int(read_max_occs_frac * n_seq_kmers))
            else:
                read_min_occs_seq, read_max_occs_seq = (read_min_occs, read_max_occs)

            if read_min_occs_seq <= seq_occs <= read_max_occs_seq:
                seqs_ids_out.add(rec.id)
                mate_cnt[rec.id[-2:] if rec.id[-2:] in ('/1', '/2') else '/0'] += 1


        _log.debug('kmer occs histogram: %s', sorted(seq_occs_hist.items()))
        _log.debug('filter_reads: %d of %d passed; mate_cnt=%s', len(seqs_ids_out), len(in_recs),
                   mate_cnt)

        return seqs_ids_out

    def binary_op(self, op, kmer_counts_1, kmer_counts_2, result_counter_cap=255):
        assert isinstance(kmer_counts_1, collections.Counter)
        assert isinstance(kmer_counts_2, collections.Counter)
        if op == 'intersect':
            result = kmer_counts_1 & kmer_counts_2
        elif op == 'union':
            result = collections.Counter({k: (kmer_counts_1[k] + kmer_counts_2[k])
                                          for k in (set(kmer_counts_1.keys()) | set(kmer_counts_2.keys()))})
        elif op == 'kmers_subtract':
            result = collections.Counter(util_misc.subdict(kmer_counts_1,
                                                           set(kmer_counts_1.keys()) - set(kmer_counts_2.keys())))
        elif op == 'counters_subtract':
            result = kmer_counts_1 - kmer_counts_2
        else:
            raise ValueError('Unknown operation: {}'.format(op))
        return self._filter_kmer_counts(counts=result, counter_cap=result_counter_cap)

# end: class KmcPy    

kmcpy = KmcPy()

def _inp(fname):
    """Return full path to a test input file for this module"""
    return os.path.join(util_file.get_test_input_path(), 'TestKmers', fname)

def _stringify(arg):
    """Return a string based on `arg`, suitable for use as a pytest test id"""
    return util_file.string_to_file_name(str(arg))

def _do_build_kmer_db(t_dir, val_cache, seq_files, kmer_db_opts):
    """Build a database of kmers from given sequence file(s) using given options.

    Args:
       t_dir: dir where to build the kmer dbase
       val_cache: a dict where results can be cached
       seq_files: a string of sequence file name(s), space-separated
       kmer_db_opts: options to build_kmer_db command, as one string
    Returns:
       an argparse.Namespace() with the following attrs:
          kmer_db: path to the kmer db created from seq_files
          kmc_kmer_counts: map from kmer to count, as computed by kmc
          kmc_db_args: the result of parsing kmer_db_opts
    """
    key = (seq_files, kmer_db_opts)
    if key in val_cache:
        _log.debug('reusing cached kmer db: %s', key)
        return val_cache[key]

    _log.debug('constructing kmer db: %s', key)

    k_db = os.path.join(t_dir, 'bld_kmer_db_{}'.format(hash(key)))
    assert not kmc.KmcTool().is_kmer_db(k_db)
    seq_files = list(map(_inp, seq_files.split()))
    kmer_db_args = util_cmd.run_cmd(module=kmer_utils, cmd='build_kmer_db',
                                    args=seq_files + [k_db] + kmer_db_opts.split() + ['--memLimitGb', 4]).args_parsed
    assert kmc.KmcTool().is_kmer_db(k_db)
    kmc_kmer_counts=kmc.KmcTool().get_kmer_counts(k_db, threads=kmer_db_args.threads)
    _log.debug('KMER_DB_FIXTURE: param=%s counts=%d db=%s', key, len(kmc_kmer_counts), k_db)
    return val_cache.setdefault(key, argparse.Namespace(kmer_db=k_db,
                                                        kmc_kmer_counts=kmc_kmer_counts,
                                                        kmer_db_args=kmer_db_args))

KMER_DBS_EMPTY = [('empty.fasta', '')]

KMER_DBS_SMALL = [
    ('ebola.fasta.gz', '-k 7'),
    ('almost-empty-2.bam', '-k 23 --singleStrand'),
    ('almost-empty-2.bam', '-k 5 --minOccs 1 --maxOccs 5 --counterCap 3'),
    ('tcgaattt.fasta', '-k 7'),
    ('ambig_bases.fasta', '-k 7'),
    ('palindromic_kmers.fasta', '-k 7')
]

KMER_DBS_MEDIUM = [
    ('test-reads.bam test-reads-human.bam', '-k 17'),
]

@pytest.fixture(scope='module')
def dict_module():
    return dict()

@pytest.fixture(scope='module', params=KMER_DBS_EMPTY+KMER_DBS_SMALL, ids=_stringify)
def kmer_db_fixture(request, tmpdir_module, dict_module):
    yield _do_build_kmer_db(tmpdir_module, dict_module, *request.param)

@pytest.fixture(scope='module', params=KMER_DBS_EMPTY+KMER_DBS_SMALL, ids=_stringify)
def kmer_db_fixture2(request, tmpdir_module, dict_module):
    yield _do_build_kmer_db(tmpdir_module, dict_module, *request.param)

@pytest.mark.parametrize("kmer_db_fixture", KMER_DBS_EMPTY+KMER_DBS_SMALL+KMER_DBS_MEDIUM,
                         ids=_stringify, indirect=["kmer_db_fixture"])
def test_build_kmer_db(kmer_db_fixture):
    _test_build_kmer_db(kmer_db_fixture)

def _test_build_kmer_db(kmer_db_fixture):
    assert kmc.KmcTool().is_kmer_db(kmer_db_fixture.kmer_db)

    kmer_db_info = kmc.KmcTool().get_kmer_db_info(kmer_db_fixture.kmer_db)
    assert kmer_db_info.kmer_size == kmer_db_fixture.kmer_db_args.kmer_size
    assert kmer_db_info.min_occs == kmer_db_fixture.kmer_db_args.min_occs
    assert kmer_db_info.max_occs == kmer_db_fixture.kmer_db_args.max_occs

    kmcpy_kmer_counts = kmcpy.compute_kmer_counts(**vars(kmer_db_fixture.kmer_db_args))
    assert kmer_db_info.total_kmers == len(kmcpy_kmer_counts)
    assert kmer_db_fixture.kmc_kmer_counts == kmcpy_kmer_counts

###########
SEQ_FILES = [ 
    # 'empty.fasta',
    # 'ebola.fasta.gz',
    # 'almost-empty-2.bam',
    # 'test-reads.bam',
    # 'test-reads-human.bam',
    # 'tcgaattt.fasta',
    'G5012.3.fasta',
    'G5012.3.mini.bam',
]

KMER_SIZES = [1, 2, 7, 17, 27, 31, 55, 63]
#KMER_SIZES = [1, 17, 31, 55]
STRAND_OPTS = ['', '--singleStrand']
KMER_OCCS_OPTS = [ '', '--minOccs 1', '--minOccs 10', '--maxOccs 10' ]
NTHREADS = [1, 2, 12]
COMBO_OPTS = [(seq_file, '-k{} {} {} --threads {}'.format(kmer_size, strand_opt, kmer_occs_opt, nthreads))
              for seq_file, kmer_size, strand_opt, kmer_occs_opt, nthreads \
              in itertools.product(SEQ_FILES, KMER_SIZES, STRAND_OPTS,
                                   KMER_OCCS_OPTS, NTHREADS)]

@pytest.mark.slow
@pytest.mark.parametrize("kmer_db_fixture", COMBO_OPTS, ids=_stringify, indirect=["kmer_db_fixture"])
def test_build_kmer_db_combo(kmer_db_fixture):
    _test_build_kmer_db(kmer_db_fixture)


##############################################################################################

def _test_filter_reads(kmer_db_fixture, reads_file, filter_opts, tmpdir_function):
    """Test read filtering.

    Args:
      kmer_size: kmer size
      single_strand: whether to canonicalize kmers
      kmers_fasta: kmers for filtering will be extracted from here
      reads_bam: reads to filter with kmers extracted from kmers_fasta

    """
    assert kmc.KmcTool().is_kmer_db(kmer_db_fixture.kmer_db)

    reads_file = _inp(reads_file)
    reads_file_out = os.path.join(tmpdir_function, 'reads_out' + util_file.uncompressed_file_type(reads_file))

    filter_args = util_cmd.run_cmd(module=kmer_utils, cmd='filter_reads',
                                   args=[kmer_db_fixture.kmer_db, 
                                         reads_file, reads_file_out] + filter_opts.split()).args_parsed

    _log.debug('Running filte: kmer_db_args=%s filter_arg=%s', kmer_db_fixture.kmer_db_args, filter_args)
    filtered_ids_expected = kmcpy.filter_reads(db_kmer_counts=kmer_db_fixture.kmc_kmer_counts,
                                               kmer_size=kmer_db_fixture.kmer_db_args.kmer_size,
                                               single_strand=kmer_db_fixture.kmer_db_args.single_strand,
                                               **vars(filter_args))

    reads_file_out_ids_txt = reads_file_out+'.ids.txt'
    read_utils.read_names(reads_file_out, reads_file_out_ids_txt)
    reads_out_ids = util_file.slurp_file(reads_file_out_ids_txt).strip().split()

    _log.debug('FILT %d %d %s %s %s', len(_list_seq_recs(reads_file)), len(_list_seq_recs(reads_file_out)),
               kmer_db_fixture.kmer_db, reads_file, filter_opts)
    def normed_read_ids(ids): return set(map(_strip_mate_num, ids))

    assert normed_read_ids(reads_out_ids) == normed_read_ids(filtered_ids_expected)

# end: def _test_filter_reads(kmer_db_fixture, reads_file, filter_opts, tmpdir_function)

@pytest.mark.parametrize("kmer_db_fixture", [('empty.fasta', '')], ids=_stringify, indirect=["kmer_db_fixture"])
@pytest.mark.parametrize("reads_file", ['empty.fasta', 'tcgaattt.fasta', 'G5012.3.subset.bam'])
@pytest.mark.parametrize("filter_opts", ['', '--readMinOccs 1', '--readMaxOccs 2'])
def test_filter_with_empty_db(kmer_db_fixture, reads_file, filter_opts, tmpdir_function):
    _test_filter_reads(**locals())

@pytest.mark.parametrize("kmer_db_fixture", [('ebola.fasta.gz', '-k 7')], ids=_stringify, indirect=["kmer_db_fixture"])
@pytest.mark.parametrize("reads_file", [pytest.param('G5012.3.testreads.bam', marks=pytest.mark.slow),
                                        'G5012.3.subset.bam'])
@pytest.mark.parametrize("filter_opts", ['--dbMinOccs 7  --readMinOccs 93',
                                         '--dbMinOccs 4 --readMinOccsFrac .6',
                                         '--readMinOccsFrac .4 --readMaxOccsFrac .55'])
def test_filter_reads(kmer_db_fixture, reads_file, filter_opts, tmpdir_function):
    _test_filter_reads(**locals())


@pytest.mark.parametrize("kmer_db_fixture", KMER_DBS_EMPTY+KMER_DBS_SMALL[:1], ids=_stringify, indirect=["kmer_db_fixture"])
@pytest.mark.parametrize("set_to_val", [1, util_misc.MAX_INT32])
def test_kmer_set_counts(kmer_db_fixture, tmpdir_function, set_to_val):
    db_with_set_counts = os.path.join(tmpdir_function, 'set_counts_db')
    util_cmd.run_cmd(module=kmer_utils, cmd='kmers_set_counts',
                     args=[kmer_db_fixture.kmer_db, set_to_val, db_with_set_counts])
    new_counts = kmc.KmcTool().get_kmer_counts(db_with_set_counts)
    assert set(new_counts.keys()) == set(kmer_db_fixture.kmc_kmer_counts.keys())
    assert not new_counts  or  set(new_counts.values()) == set([set_to_val])


@pytest.mark.parametrize("op", ('intersect', 'union', 'kmers_subtract', 'counters_subtract'))
def test_kmers_binary_op(kmer_db_fixture, kmer_db_fixture2, op, tmpdir_function):
    if kmer_db_fixture.kmer_db_args.kmer_size != kmer_db_fixture2.kmer_db_args.kmer_size:
        pytest.skip('(always skip) binary ops not defined on kmers of different size')
    db_result = os.path.join(tmpdir_function, 'op_result')
    _log.debug('fixture1args=%s', kmer_db_fixture.kmer_db_args)
    _log.debug('fixture2args=%s', kmer_db_fixture2.kmer_db_args)
    _log.debug('op=%s', op)
    args = util_cmd.run_cmd(module=kmer_utils, cmd='kmers_binary_op',
                            args=[op, kmer_db_fixture.kmer_db, kmer_db_fixture2.kmer_db, db_result]).args_parsed

    kmc_counts = kmc.KmcTool().get_kmer_counts(db_result)
    kmcpy_counts = kmcpy.binary_op(op, kmer_db_fixture.kmc_kmer_counts, kmer_db_fixture2.kmc_kmer_counts,
                                   result_counter_cap=args.result_counter_cap)

    assert kmc_counts == kmcpy_counts
