"""Unit tests for kmers.py"""

__author__ = "ilya@broadinstitute.org"

import os
import collections
import argparse
import inspect

from test import assert_equal_contents, assert_equal_bam_reads, make_slow_test_marker

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pytest

import kmers
import util.cmd
import util.file
import util.misc
import tools.kmc
import tools.samtools

slow_test = make_slow_test_marker()  # pylint: disable=invalid-name

#################################
# Some general utils used below #
#################################

def locals_sans_self():
    return {k: v for k, v in inspect.currentframe().f_back.f_locals.items() if k != 'self'}

def locals_dict(vars):
    return {k: v for k, v in inspect.currentframe().f_back.f_locals.items() if k in vars.split()}

def _seq_as_str(s):  # pylint: disable=invalid-name
    """Return a sequence as a str, regardless of whether it was a str, a Seq or a SeqRecord"""
    if isinstance(s, Seq):
        return str(s)
    if isinstance(s, SeqRecord):
        return str(s.seq)
    return s

def _yield_seqs_as_strs(seqs):
    """Yield sequence(s) from `seqs` as strs.  seqs can be a str/SeqRecord/Seq, a filename of a sequence file,
    or an iterable of these.  If a filename of a sequence file, all sequences from that file are yielded."""
    for seq in util.misc.make_seq(seqs, (str, SeqRecord, Seq)):
        seq = _seq_as_str(seq)
        if not any(seq.endswith(ext) for ext in '.fasta .fasta.gz .fastq .fastq.gz .bam'):
            yield seq
        else:
            with util.file.tmp_dir(suffix='_seqs_as_strs') as t_dir:
                if seq.endswith('.bam'):
                    t_fa = os.path.join(t_dir, 'bam2fa.fasta')
                    tools.samtools.SamtoolsTool().bam2fa(seq, t_fa)
                    seq = t_fa
                with util.file.open_or_gzopen(seq, 'rt') as seq_f:
                    for rec in Bio.SeqIO.parse(seq_f, util.file.uncompressed_file_type(seq)[1:]):
                        yield str(rec.seq)

def _getargs(args, valid_args):
    """Extract valid args from an argparse.Namespace or a dict.  Returns a dict containing keys from `args`
    that are in `valid_args`; `valid_args` is a space-separated list of valid args."""
    return util.misc.subdict(vars(args) if isinstance(args, argparse.Namespace) else args, valid_args.split())

#################################################################################################################

class KmcPy(object):
    """Reimplementation of some kmc functions in simple Python code."""

    #
    # To help generate the expected correct output, we reimplement in simple Python code some
    # of KMC's functionality.   This is also to make up for KMC's lack of a public test suite
    # ( https://github.com/refresh-bio/KMC/issues/55 ).
    #

    def _revcomp(self, kmer):
        """Return the reverse complement of a kmer, given as a string"""
        return str(Seq(kmer, IUPAC.unambiguous_dna).reverse_complement())

    def _canonicalize(self, kmer):
        """Return the canonical version of a kmer"""
        return min(kmer, self._revcomp(kmer))

    def _compute_kmers(self, seq_strs, kmer_size, single_strand, **ignore):
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
                if seq[i] not in 'TCGA':
                    invalid_lo = max(0, i-kmer_size+1)
                    invalid_hi = min(n_kmers, i+1)
                    valid_kmer[invalid_lo:invalid_hi] = [False]*(invalid_hi-invalid_lo)

            for i in range(n_kmers):
                if valid_kmer[i]:
                    kmer = seq[i:i+kmer_size]
                    yield kmer if single_strand else self._canonicalize(kmer)

    def compute_kmer_counts(self, seq_files, kmer_size, min_occs=None, max_occs=None,
                            counter_cap=None, single_strand=False, **ignore):
        """Yield kmer counts of seq(s).  Unless `single_strand` is True, each kmer is
        canonicalized before being counted.  Kmers containing non-TCGA bases are skipped.
        Kmers with fewer than `min_occs` or more than `max_occs` occurrences
        are dropped, and kmer counts capped at `counter_cap`, if these args are given.
        """
        counts = collections.Counter(self._compute_kmers(_yield_seqs_as_strs(seq_files), kmer_size, single_strand))
        if any((min_occs, max_occs, counter_cap)):
            counts = dict((kmer, min(count, counter_cap or count)) \
                          for kmer, count in counts.items() \
                          if (count >= (min_occs or count)) and \
                          (count <= (max_occs or count)))
        return counts

    def filter_seqs(self, db_kmer_counts, seq_files, kmer_size, single_strand,
                    read_min_occs=None, read_max_occs=None, **ignore):
        seqs_out = []
        read_min_occs, read_max_occs = tools.kmc.KmcTool()._infer_filter_reads_params(read_min_occs, read_max_occs) # pylint: disable=protected-access

        for seq in _yield_seqs_as_strs(seq_files):
            seq_kmer_counts = self.compute_kmer_counts(seq, kmer_size=kmer_size, single_strand=single_strand)
            seq_occs = sum([seq_count for kmer, seq_count in seq_kmer_counts.items() \
                            if kmer in db_kmer_counts])
            read_min_occs_seq, read_max_occs_seq = map(lambda v: int(v*(len(seq)-kmer_size+1)) # pylint: disable=cell-var-from-loop
                                                       if isinstance(v, float) else v,
                                                       (read_min_occs, read_max_occs))

            if (read_min_occs_seq or seq_occs) <= seq_occs <= (read_max_occs_seq or seq_occs):
                seqs_out.append(seq)

        return seqs_out

# end: class KmcPy    

kmcpy = KmcPy()

def _inp(fname):
    """Return full path to a test input file for this module"""
    return os.path.join(util.file.get_test_input_path(), 'TestKmers', fname)

def _stringify(par): 
    return util.file.string_to_file_name(str(par))

BUILD_KMER_DB_TESTS = [
    ('empty.fasta', ''),
    ('ebola.fasta.gz', '-k 7 --threads 2'),
    ('almost-empty-2.bam', '-k 23 --singleStrand'),
    ('almost-empty-2.bam', '-k 5 --minOccs 1 --maxOccs 5 --counterCap 3 --threads 2'),
    ('test-reads.bam test-reads-human.bam', '-k 17'),
    pytest.param(('tcgaattt', ' -k 7 --threads 11'),
                 marks=(pytest.mark.skipif(util.misc.available_cpu_count() < 11,
                                           reason='needs 11+ threads to show bug'),
                        pytest.mark.xfail(reason='kmc bug')))
]

@pytest.fixture(scope='module',
                params=BUILD_KMER_DB_TESTS,
                ids=_stringify)
def kmer_db_fixture(request, tmpdir_module):
    """Build a database of kmers from given sequence file(s)"""
    k_db = os.path.join(tmpdir_module, 'bld_kmer_db')
    seq_files, opts = request.param
    seq_files = list(map(_inp, seq_files.split()))
    kmer_db_args = util.cmd.run_cmd(module=kmers, cmd='build_kmer_db',
                                    args=seq_files + [k_db] + opts.split() + ['--memLimitGb', 4]).args_parsed
    kmc_kmer_counts=tools.kmc.KmcTool().get_kmer_counts(k_db, threads=kmer_db_args.threads)

    yield argparse.Namespace(kmer_db=k_db,
                             kmc_kmer_counts=kmc_kmer_counts,
                             kmer_db_args = kmer_db_args)

def test_build_kmer_db(kmer_db_fixture):
    assert tools.kmc.KmcTool().is_kmer_db(kmer_db_fixture.kmer_db)

    kmer_db_info = tools.kmc.KmcTool().get_kmer_db_info(kmer_db_fixture.kmer_db)
    assert kmer_db_info.kmer_size == kmer_db_fixture.kmer_db_args.kmer_size
    assert kmer_db_info.min_occs == kmer_db_fixture.kmer_db_args.min_occs
    assert kmer_db_info.max_occs == kmer_db_fixture.kmer_db_args.max_occs

    kmcpy_kmer_counts = kmcpy.compute_kmer_counts(**vars(kmer_db_fixture.kmer_db_args))
    assert kmer_db_info.total_kmers == len(kmcpy_kmer_counts)
    assert kmer_db_fixture.kmc_kmer_counts == kmcpy_kmer_counts

def _test_filter_by_kmers(kmer_db_fixture, reads_file, filter_opts, tmpdir_function):
    """Test read filtering.

    Args:
      kmer_size: kmer size
      single_strand: whether to canonicalize kmers
      kmers_fasta: kmers for filtering will be extracted from here
      reads_bam: reads to filter with kmers extracted from kmers_fasta

    """
    assert tools.kmc.KmcTool().is_kmer_db(kmer_db_fixture.kmer_db)

    reads_file = _inp(reads_file)
    reads_file_out = os.path.join(tmpdir_function, 'reads_out' + util.file.uncompressed_file_type(reads_file))

    filter_args = util.cmd.run_cmd(module=kmers, cmd='filter_by_kmers',
                                   args=[kmer_db_fixture.kmer_db, 
                                         reads_file, reads_file_out] + filter_opts.split()).args_parsed


    filtered_seqs_expected = kmcpy.filter_seqs(db_kmer_counts=kmer_db_fixture.kmc_kmer_counts, 
                                               seq_files=reads_file,
                                               kmer_size=kmer_db_fixture.kmer_db_args.kmer_size,
                                               single_strand=kmer_db_fixture.kmer_db_args.single_strand,
                                               **vars(filter_args))
    assert set(_yield_seqs_as_strs(reads_file_out)) == set(filtered_seqs_expected)

# end: def _test_filter_by_kmers(kmer_db_fixture, reads_file, filter_opts, tmpdir_function)

@pytest.mark.parametrize("reads_file,filter_opts",
                         [('ebola.fasta', '')], ids=_stringify)
def test_filter_by_kmers(kmer_db_fixture, reads_file, filter_opts, tmpdir_function):
    _test_filter_by_kmers(**locals())
