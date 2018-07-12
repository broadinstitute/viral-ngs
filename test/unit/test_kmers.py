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
            with util.file.tmp_dir(suffix='seqs_as_strs') as t_dir:
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

    @staticmethod
    def _revcomp(kmer):
        """Return the reverse complement of a kmer, given as a string"""
        return str(Seq(kmer, IUPAC.unambiguous_dna).reverse_complement())

    @staticmethod
    def _canonicalize(kmer):
        """Return the canonical version of a kmer"""
        return min(kmer, KmcPy._revcomp(kmer))

    @staticmethod
    def _compute_kmers(seq_strs, kmer_size, single_strand, **ignore):
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
                    yield kmer if single_strand else KmcPy._canonicalize(kmer)

    @staticmethod
    def compute_kmer_counts(seq_files, kmer_size, min_occs=None, max_occs=None,
                            counter_cap=None, single_strand=False, **kwargs):
        """Yield kmer counts of seq(s).  Unless `single_strand` is True, each kmer is
        canonicalized before being counted.  Kmers containing non-TCGA bases are skipped.
        Kmers with fewer than `min_occs` or more than `max_occs` occurrences
        are dropped, and kmer counts capped at `counter_cap`, if these args are given.
        """
        counts = collections.Counter(KmcPy._compute_kmers(_yield_seqs_as_strs(seq_files), kmer_size, single_strand))
        if any((min_occs, max_occs, counter_cap)):
            counts = dict((kmer, min(count, counter_cap or count)) \
                          for kmer, count in counts.items() \
                          if (count >= (min_occs or count)) and \
                          (count <= (max_occs or count)))
        return counts

# end: class KmcPy    

def _inp(fname):
    return os.path.join(util.file.get_test_input_path(), 'TestKmers', fname)

def _stringify(par): 
    return util.file.string_to_file_name(str(par))

BUILD_KMER_DB_TESTS = [
    ('empty.fasta', '-k 21'),
    ('ebola.fasta.gz', '-k 7 --threads 2'),
    ('almost-empty-2.bam', '-k 23 --singleStrand'),
    ('almost-empty-2.bam', '-k 5 --minOccs 1 --maxOccs 5 --counterCap 3 --threads 2'),
    ('test-reads.bam test-reads-human.bam', '-k 17'),
    pytest.param(('tcgaattt', ' -k 7 --threads 11'),
                 marks=(pytest.mark.skipif(util.misc.available_cpu_count() < 11,
                                           reason='needs 11+ threads to show bug'),
                        pytest.mark.xfail(reason='kmc bug 1')))
]


@pytest.fixture(scope='module',
                params=BUILD_KMER_DB_TESTS,
                ids=_stringify)
def kmer_db_fixture(request):
    """Build a database of kmers from given sequence file(s)"""
    with util.file.tmp_dir(suffix='build_kmer_db') as t_dir:
        k_db = os.path.join(t_dir, 'bld_kmer_db')
        seq_files, opts = request.param
        seq_files = list(map(_inp, seq_files.split()))
        kmer_db_args = util.cmd.run_cmd(module=kmers, cmd='build_kmer_db', 
                                        args=seq_files + [k_db] + opts.split() + ['--memLimitGb', 4]).args_parsed

        kmers_txt = k_db+'.kmer_counts.txt'
        kmers.dump_kmer_counts(kmer_db=k_db, out_kmers=kmers_txt, threads=kmer_db_args.threads)

        return argparse.Namespace(kmer_db=k_db, kmc_kmer_counts=tools.kmc.KmcTool().read_kmer_counts(kmers_txt),
                                  kmcpy_kmer_counts=KmcPy.compute_kmer_counts(**vars(kmer_db_args)))

def test_build_kmer_db(kmer_db_fixture):
    assert kmer_db_fixture.kmc_kmer_counts == kmer_db_fixture.kmcpy_kmer_counts


class TestKmers(object):

    """Test commands for manipulating kmers"""

    # to test:
    #   empty bam, empty fasta, empty both
    #   getting kmers from bam (here just convert the fasta to it?), from fastq,
    #   from multiple files.
    #   ambiguity codes, gaps, Ns
    #   bams with different mixes of read groups, mix of paired/unpaired etc

    def _filter_seqs(self, db_kmer_counts, seqs, kmer_size, single_strand, read_min_occs=None, read_max_occs=None):
        seqs_out = []
        read_min_occs, read_max_occs = tools.kmc.KmcTool()._infer_filter_reads_params(read_min_occs, read_max_occs) # pylint: disable=protected-access

        for seq in util.misc.make_seq(seqs, (str, SeqRecord, Seq)):
            seq_kmer_counts = KmcPy.compute_kmer_counts(seq, kmer_size=kmer_size, single_strand=single_strand)
            seq_occs = sum([seq_count for kmer, seq_count in seq_kmer_counts.items() \
                            if kmer in db_kmer_counts])
            read_min_occs_seq, read_max_occs_seq = map(lambda v: int(v*(len(seq)-kmer_size+1)) # pylint: disable=cell-var-from-loop
                                                       if isinstance(v, float) else v,
                                                       (read_min_occs, read_max_occs))
            #print('seq=', seq, 'kmer_counts=', seq_kmer_counts, 'seq_occs=', seq_occs)

            if (read_min_occs_seq or seq_occs) <= seq_occs <= (read_max_occs_seq or seq_occs):
                seqs_out.append(seq)

        return seqs_out

    @pytest.mark.parametrize("kmer_size", [17, 31])
    @pytest.mark.parametrize("single_strand", [False, True])
    @pytest.mark.parametrize("read_min_occs", [80, .7])
    @pytest.mark.parametrize("kmers_fasta,reads_bam", [
        ('ebola.fasta', 'G5012.3.subset.bam'),
    ])
    def test_read_filtering(self, kmer_size, single_strand, kmers_fasta, reads_bam, read_min_occs):
        self._test_read_filtering(kmer_size, single_strand, kmers_fasta, reads_bam, read_min_occs)

    @pytest.mark.xfail(reason="kmc bug, reported at https://github.com/refresh-bio/KMC/issues/86")
    @pytest.mark.parametrize("args", [dict(kmer_size=1, single_strand=False, read_min_occs=1, kmers_fasta='empty.fasta',
                                           reads_bam='empty.bam')])
    def test_read_filtering_fail_1(self, args):
        self._test_read_filtering(**args)

    def _test_read_filtering(self, kmer_size, single_strand, kmers_fasta, reads_bam,
                             read_min_occs=None, read_max_occs=None):
        """Test read filtering.

        Args:
          kmer_size: kmer size
          single_strand: whether to canonicalize kmers
          kmers_fasta: kmers for filtering will be extracted from here
          reads_bam: reads to filter with kmers extracted from kmers_fasta

        """
        with util.file.tmp_dir(suffix='kmctest_reafilt') as t_dir:
            def _cmdflag(flag, val):
                return [] if val in (None, False) else ([flag] if val is True else [flag, val])

            kmers_fasta = os.path.join(util.file.get_test_input_path(), kmers_fasta)
            kmer_db = os.path.join(t_dir, 'kmer_db')
            kmer_db_args = util.cmd.run_cmd(kmers, 'build_kmer_db',
                                            ['--memLimitGb', 4, '-k', kmer_size, kmers_fasta, kmer_db]
                                            +_cmdflag('--singleStrand', single_strand)).args_parsed

            kmers_fasta_seqs = list(Bio.SeqIO.parse(kmers_fasta, 'fasta'))
            kmers_fasta_kmer_counts = KmcPy.compute_kmer_counts(kmers_fasta_seqs,
                                                                        **_getargs(kmer_db_args,
                                                                                   'kmer_size single_strand '
                                                                                   'min_occs max_occs counter_cap'))

            kmers_txt = os.path.join(t_dir, 'kmers.txt')
            util.cmd.run_cmd(kmers, 'dump_kmer_counts', [kmer_db, kmers_txt])

            kmc_kmer_counts = tools.kmc.KmcTool().read_kmer_counts(kmers_txt)

            assert len(kmc_kmer_counts) == len(kmers_fasta_kmer_counts)
            assert kmc_kmer_counts == kmers_fasta_kmer_counts

            reads_bam = os.path.join(util.file.get_test_input_path(), reads_bam)
            with tools.samtools.SamtoolsTool().bam2fa_tmp(reads_bam) as (reads_1, reads_2):
                reads_1_filt = os.path.join(t_dir, 'ebola.reads.1.filt.fasta')
                util.cmd.run_cmd(kmers, 'filter_by_kmers', [kmer_db, reads_1, reads_1_filt] \
                                                            + _cmdflag('--readMinOccs', read_min_occs) \
                                                            + _cmdflag('--readMaxOccs', read_max_occs))
                reads_1_seqs = tuple(Bio.SeqIO.parse(reads_1, 'fasta'))

                reads_1_filt_expected = self._filter_seqs(kmc_kmer_counts, reads_1_seqs,
                                                          kmer_size=kmer_size, single_strand=single_strand,
                                                          read_min_occs=read_min_occs)
                reads_1_filt_seqs = tuple(Bio.SeqIO.parse(reads_1_filt, 'fasta'))
                #assert 0 < len(reads_1_filt_seqs) < len(reads_1_seqs)
                def SeqRecord_data(r):  # pylint: disable=invalid-name
                    return (r.id, r.seq)
                assert  sorted(map(SeqRecord_data, reads_1_filt_seqs)) == \
                    sorted(map(SeqRecord_data, reads_1_filt_expected))

    # end: def _test_read_filtering(self, kmer_size, single_strand, kmers_fasta, reads_bam,
    #                               read_min_occs=None, read_max_occs=None)

    def _test_filter_by_kmers(kmer_db, in_reads, out_reads, db_min_occs=None, db_max_occs=None,
                              read_min_occs=None, read_max_occs=None, hard_mask=False, threads=None):
        """Test the operation of filter_by_kmers"""
        pass


# end: class TestKmers(object)
