# Unit tests for intrahost.py

__author__ = "dpark@broadinstitute.org"

# built-ins
from collections import OrderedDict
import os
import os.path
import shutil
import tempfile
import itertools
import argparse
import unittest

# third-party
import Bio
import Bio.SeqRecord
import Bio.Seq

# module-specific
import intrahost
import util.file
import util.vcf
import test
from intrahost import AlleleFieldParser
import interhost
import tools.mafft


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in intrahost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


def makeTempFasta(seqs):
    fn = util.file.mkstempfname('.fasta')
    with open(fn, 'wt') as outf:
        for line in util.file.fastaMaker(seqs):
            outf.write(line)
    return fn


class MockVphaserOutput:
    ''' This creates test data that pretends to be the output from
        tools.vphaser2.Vphaser2Tool.iterate
    '''

    def __init__(self):
        self.isnvs = {}
        self.chroms = []

    def add_snp(self, chrom, pos, acounts, libinfo=None):
        # Add an iSNP at this chrom,pos. acounts is a list of triples:
        # (allele, fwd count, rev count).
        assert isinstance(pos, int) and pos > 0 and len(acounts) > 1
        for a, f, r in acounts:
            assert a in ('A', 'C', 'G', 'T')
            assert f >= 0 and r >= 0 and f + r > 0

        # Add libinfo to acounts
        if libinfo is None:
            acounts = [[a, f, r, 1, [[f, r]]] for a, f, r in acounts]
        else:
            # ...
            raise Exception

        acounts = reversed(sorted((acount[1] + acount[2],) + tuple(acount) for acount in acounts))
        acounts = [tuple(acount[1:]) for acount in acounts]
        if chrom not in self.chroms:
            self.chroms.append(chrom)
        self.isnvs.setdefault(chrom, {})
        self.isnvs[chrom].setdefault(pos, {})
        self.isnvs[chrom][pos]['snp'] = acounts

    def add_indel(self, chrom, pos, acounts, libinfo=None):
        ''' Add an iSNP at this chrom,pos. acounts is a list of triples:
            (allele, fwd count, rev count). allele is simply a small
            sequence with no "I" or "D" prefixing it (we'll add the I/D).
        '''
        assert isinstance(pos, int) and pos > 0 and len(acounts) > 1
        for a, f, r in acounts:
            assert isinstance(a, str)
            assert f >= 0 and r >= 0 and f + r > 0
        assert '' in set(a for a, f, r in acounts)

        # Add libinfo to acounts
        if libinfo is None:
            acounts = [[a, f, r, 1, [[f, r]]] for a, f, r in acounts]
        else:
            # ...
            raise Exception

        acounts = reversed(sorted((acount[1] + acount[2],) + tuple(acount) for acount in acounts))
        acounts = [list(acount[1:]) for acount in acounts]

        # vphaser funniness here
        if acounts[0][0] == '':
            # this is a set of insertions against the consensus
            acounts[0][0] = 'd'
            for i in range(1, len(acounts)):
                acounts[i][0] = 'I' + acounts[i][0]
        else:
            # this is a deletion against the consensus
            for i in range(1, len(acounts)):
                n_deleted = len(acounts[0][0]) - len(acounts[i][0])
                assert n_deleted > 0
                assert acounts[0][0][n_deleted:] == acounts[i][0]
                acounts[i][0] = 'D' + str(n_deleted)
            acounts[0][0] = 'i'

        if chrom not in self.chroms:
            self.chroms.append(chrom)
        self.isnvs.setdefault(chrom, {})
        self.isnvs[chrom].setdefault(pos, {})
        self.isnvs[chrom][pos]['lp'] = acounts

    def __iter__(self):
        for c in self.chroms:
            for p in sorted(self.isnvs[c].keys()):
                for model in self.isnvs[c][p].keys():
                    acounts = self.isnvs[c][p][model]
                    mac = sum(acount[1] + acount[2] for acount in acounts[1:])
                    tot = sum(acount[1] + acount[2] for acount in acounts)
                    yield [c, str(p), acounts[1][0], acounts[0][0],
                           '0.5', model, str(float(mac) / tot * 100.0)] \
                        + [str(AlleleFieldParser(None, *acount))
                           for acount in acounts]


class TestIntrahostFilters(unittest.TestCase):
    ''' This tests step 1 of the iSNV calling process
        (intrahost.vphaser_one_sample), which runs V-Phaser2 on
        a single sample, reformats the output slightly, and performs
        strand-bias filtering and adds library-bias statistics.
        These unit tests mock the vphaser tool output and just test
        the filtering/statistics/etc.
    '''

    def test_single_strand_bias_hard_filter(self):
        data = MockVphaserOutput()
        data.add_snp('c1', 100, [('A', 10, 20), ('T', 5, 2), ('C', 30, 500), ('G', 60, 40)])
        data.add_snp('c2', 100, [('C', 10, 2), ('T', 2, 8)])
        output = list(intrahost.filter_strand_bias(data))
        expected = ['c1', '100', 'A', 'G', None, 'snp', 23.076923076923078, 'G:60:40:60:40:1', 'A:10:20:10:20:1']
        self.assertEqual(len(output), 1)
        self.assertEqual(output[0][:4], expected[:4])
        self.assertEqual(output[0][5], expected[5])
        self.assertAlmostEqual(float(output[0][6]), expected[6], places=4)
        self.assertEqual(output[0][7:], expected[7:])


#@unittest.skipIf(tools.is_osx(), "vphaser2 osx binary from bioconda has issues")
class TestPerSample(test.TestCaseWithTmp):
    ''' This tests step 1 of the iSNV calling process
        (intrahost.vphaser_one_sample), which runs V-Phaser2 on
        a single sample, reformats the output slightly, and performs
        strand-bias filtering and adds library-bias statistics.
    '''

    def test_vphaser_one_sample(self):
        # Files here were created as follows:
        # - in.bam was copied from input directory for TestVPhaser2; see notes
        #   there on how it was created.
        # - ref.fasta was created by making two identical chromosomes, chr1
        #   and chr2, with the sequence from West Nile virus isolate
        #   WNV-1/US/BID-V7821/2011. That genome is a guess for the reference
        #   of the V-Phaser 2 test file because BLAST matches the reads to
        #   West Nile virus and that isolate has the size reported in the bam file.
        #   Note that ref.fasta is not exactly the consensus of the reads in in.bam;
        #   for example, pos 660 is C in ref.fasta, but more reads have T there
        #   than C in in.bam. So we are actually testing the case that
        #   V-Phaser 2 consensus != our consensus.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.bam')
        refFasta = os.path.join(myInputDir, 'ref.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=test._CPUS, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_expected.txt')
        self.assertEqualContents(outTab, expected)


class VcfMergeRunner:
    ''' This creates test data and feeds it to intrahost.merge_to_vcf
    '''

    def __init__(self, ref_genome=None):
        self.genomes = {}  # {sample : {chrom : bases, ...}, ...}
        self.genomeFastas = OrderedDict()  # {sample: fastaFileName, ...
        self.alignedFastas = []  # [chr1, chr2, ...]
        self.isnvs = {}
        self.sample_order = []
        self.sequence_order = OrderedDict()  # { sample: [seq1, seq2, ... ], ... }
        if ref_genome:
            self.set_ref(ref_genome)

    def set_ref(self, genome):
        self.ref = makeTempFasta(genome)

    def add_genome(self, sample_name, genome):
        self.genomes[sample_name] = dict(genome)
        self.genomeFastas[sample_name] = makeTempFasta(genome)
        if sample_name not in self.sample_order:
            self.sample_order.append(sample_name)
        self.isnvs.setdefault(sample_name, MockVphaserOutput())

    def add_snp(self, sample, chrom, pos, acounts, libinfo=None):
        assert sample in self.genomeFastas
        assert chrom in self.genomes[sample]
        assert 1 <= pos <= len(self.genomes[sample][chrom])
        assert self.genomes[sample][chrom][pos - 1] in [a for a, f, r in acounts]
        self.isnvs[sample].add_snp(chrom, pos, acounts, libinfo)

    def add_indel(self, sample, chrom, pos, acounts, libinfo=None):
        assert sample in self.genomeFastas
        assert sample in self.genomeFastas
        assert chrom in self.genomes[sample]
        assert 1 <= pos <= len(self.genomes[sample][chrom])
        if acounts[0][0] != '':  # deletion
            assert self.genomes[sample][chrom][pos - 1:].startswith(acounts[0][0])
        self.isnvs[sample].add_indel(chrom, pos, acounts, libinfo)

    def dump_isnv_tmp_file(self, sample):
        fn = util.file.mkstempfname('.txt')
        with open(fn, 'wt') as outf:
            for row in self.isnvs[sample]:
                outf.write('\t'.join(map(str, row)) + '\n')
        return fn

    def run_and_get_vcf_rows(self, retree=1, omit_samplenames=False):
        outVcf = util.file.mkstempfname('.vcf.gz')

        self.multi_align_samples(retree=retree)

        if not omit_samplenames:
            intrahost.merge_to_vcf(self.ref, outVcf, self.sample_order, list(self.dump_isnv_tmp_file(s) for s in self.sample_order),
                              self.alignedFastas)
        else:
            intrahost.merge_to_vcf(self.ref, outVcf, [], list(self.dump_isnv_tmp_file(s) for s in self.sample_order),
                              self.alignedFastas)


        with util.vcf.VcfReader(outVcf) as vcf:
            rows = list(vcf.get())
        return rows

    def multi_align_samples(self, retree=1):

        # store here a dict, {sample_name: [sequence_id1, sequence_id2, ...]}

        # because in Python3 OrderedDict does not have iteritems(), it simply has items()
        if hasattr(self.genomeFastas, "iteritems"):
            genomeKVIterator = self.genomeFastas.iteritems()
        else:
            genomeKVIterator = self.genomeFastas.items()

        for sampleName, fastaFile in genomeKVIterator:
            with util.file.compressed_open(fastaFile, 'rt') as inf:
                for seq in Bio.SeqIO.parse(inf, 'fasta'):
                    self.sequence_order.setdefault(sampleName, default=[])
                    self.sequence_order[sampleName].append(seq.id)

        inputFastas = []
        inputFastas.append(self.ref)
        inputFastas.extend(self.genomeFastas.values())
        transposedFiles = interhost.transposeChromosomeFiles(inputFastas)

        # since the FASTA files are
        for idx, filePath in enumerate(transposedFiles):

            outFile = util.file.mkstempfname('.fasta')
            outFilePath = os.path.dirname(outFile)

            alignedOutFile = tools.mafft.MafftTool().execute(
                inFastas=[os.path.abspath(filePath)],
                outFile=os.path.join(outFilePath, "{}{}.fasta".format("aligned", idx)),
                localpair=False,
                globalpair=True,
                preservecase=True,
                reorder=None,
                gapOpeningPenalty=None,
                offset=None,
                verbose=False,
                outputAsClustal=None,
                maxiters=1000,
                threads=test._CPUS,
                retree=retree)
            self.alignedFastas.append(alignedOutFile)


class TestVcfMerge(test.TestCaseWithTmp):
    ''' This tests step 2 of the iSNV calling process
        (intrahost.merge_to_vcf), which gets really nasty and tricky
        and has lots of edge cases. These unit tests mock the vphaser
        tool output and just tests the merge and VCF stuff.
    '''

    def test_empty_output(self):
        ref = makeTempFasta([('ref1', 'ATCGCA')])
        s1 = makeTempFasta([('s1-1', 'ATCGCA')])
        emptyfile = util.file.mkstempfname('.txt')
        outVcf = util.file.mkstempfname('.vcf')
        #intrahost.merge_to_vcf(ref, outVcf, ['s1'], [emptyfile], [s1])
        self.assertRaises(LookupError, intrahost.merge_to_vcf, ref, outVcf, ['s1'], [emptyfile], [s1])
        self.assertGreater(os.path.getsize(outVcf), 0)
        with util.file.compressed_open(outVcf, 'rt') as inf:
            for line in inf:
                self.assertTrue(line.startswith('#'))
        outVcf = util.file.mkstempfname('.vcf.gz')
        #intrahost.merge_to_vcf(ref, outVcf, ['s1'], [emptyfile], [s1])
        self.assertRaises(LookupError, intrahost.merge_to_vcf, ref, outVcf, ['s1'], [emptyfile], [s1])
        self.assertGreater(os.path.getsize(outVcf), 0)
        with util.file.compressed_open(outVcf, 'rt') as inf:
            for line in inf:
                self.assertTrue(line.startswith('#'))

    def test_headers_with_two_samps(self):
        ref = makeTempFasta([('ref1', 'ATCGTTCA'), ('ref2', 'GGCCC')])
        s1 = makeTempFasta([('s1-1', 'ATCGCA'), ('s1-2', 'GGCCC')])
        s2 = makeTempFasta([('s2-1', 'ATCGTTCA'), ('s2-2', 'GGCCC')])
        emptyfile = util.file.mkstempfname('.txt')
        outVcf = util.file.mkstempfname('.vcf.gz')
        self.assertRaises(
            LookupError, intrahost.merge_to_vcf, ref, outVcf, [
                's1', 's2'
            ], [
                emptyfile, emptyfile
            ], [
                s1, s2
            ])
        with util.vcf.VcfReader(outVcf) as vcf:
            self.assertEqual(vcf.samples(), ['s1', 's2'])
            self.assertEqual(vcf.chrlens(), {'ref1': 8, 'ref2': 5})

    def test_simple_snps(self):
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGGAC')])
        # ATCGGAC-
        merger.add_genome('s2', [('s2-1',  'TCGGACT')])
        # -TCGGACT
        merger.add_genome('s3', [('s3-1',  'TCGGACT')])
        # -TCGGACT
        merger.add_snp('s1', 's1-1', 3, [('C', 80, 80), ('A', 20, 20)])
        merger.add_snp('s2', 's2-1', 2, [('C', 90, 90), ('A', 10, 10)])
        merger.add_snp('s3', 's3-1', 5, [('A', 70, 70), ('T', 30, 30)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 3)
        self.assertEqual(rows[0].ref, 'C')
        self.assertEqual(rows[0].alt, 'A')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.2')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '0:0.1')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '0:0.0')
        self.assertEqual(rows[1].contig, 'ref1')
        self.assertEqual(rows[1].pos + 1, 6)
        self.assertEqual(rows[1].ref, 'A')
        self.assertEqual(rows[1].alt, 'T')
        self.assertEqual(':'.join(rows[1][0].split(':')[:2]), '0:0.0')
        self.assertEqual(':'.join(rows[1][1].split(':')[:2]), '0:0.0')
        self.assertEqual(':'.join(rows[1][2].split(':')[:2]), '0:0.3')

    def test_simple_snps_guess_samplenames(self):
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGGAC')])
        # ATCGGAC-
        merger.add_genome('s2', [('s2-1',  'TCGGACT')])
        # -TCGGACT
        merger.add_genome('s3', [('s3-1',  'TCGGACT')])
        # -TCGGACT
        merger.add_snp('s1', 's1-1', 3, [('C', 80, 80), ('A', 20, 20)])
        merger.add_snp('s2', 's2-1', 2, [('C', 90, 90), ('A', 10, 10)])
        merger.add_snp('s3', 's3-1', 5, [('A', 70, 70), ('T', 30, 30)])
        rows = merger.run_and_get_vcf_rows(omit_samplenames=True)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 3)
        self.assertEqual(rows[0].ref, 'C')
        self.assertEqual(rows[0].alt, 'A')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.2')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '0:0.1')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '0:0.0')
        self.assertEqual(rows[1].contig, 'ref1')
        self.assertEqual(rows[1].pos + 1, 6)
        self.assertEqual(rows[1].ref, 'A')
        self.assertEqual(rows[1].alt, 'T')
        self.assertEqual(':'.join(rows[1][0].split(':')[:2]), '0:0.0')
        self.assertEqual(':'.join(rows[1][1].split(':')[:2]), '0:0.0')
        self.assertEqual(':'.join(rows[1][2].split(':')[:2]), '0:0.3')

    def test_snps_with_varying_read_depth(self):
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGGAC')])
        # ATCGGAC-
        merger.add_genome('s2', [('s2-1', 'TCGGACT')])
        # -TCGGACT
        merger.add_genome('s3', [('s3-1', 'TCGGACT')])
        # -TCGGACT
        merger.add_snp('s1', 's1-1', 3, [('C', 60, 100), ('A', 25, 15)])
        merger.add_snp('s2', 's2-1', 2, [('C', 12, 6), ('A', 2, 0)])
        merger.add_snp('s3', 's3-1', 5, [('A', 7, 4), ('T', 2, 3)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 3)
        self.assertEqual(rows[0].ref, 'C')
        self.assertEqual(rows[0].alt, 'A')
        self.assertEqual(':'.join(rows[0][0].split(':')[:3]), '0:0.2:200')
        self.assertEqual(':'.join(rows[0][1].split(':')[:3]), '0:0.1:20')
        self.assertEqual(':'.join(rows[0][2].split(':')[:3]), '0:0.0:.')
        self.assertEqual(rows[1].contig, 'ref1')
        self.assertEqual(rows[1].pos + 1, 6)
        self.assertEqual(rows[1].ref, 'A')
        self.assertEqual(rows[1].alt, 'T')
        self.assertEqual(':'.join(rows[1][0].split(':')[:3]), '0:0.0:.')
        self.assertEqual(':'.join(rows[1][1].split(':')[:3]), '0:0.0:.')
        self.assertEqual(':'.join(rows[1][2].split(':')[:3]), '0:0.3125:16')

    def test_snps_downstream_of_indels(self):
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGTTGACT')])
        merger.add_genome('s2', [('s2-1', 'TCGTTGACT')])
        merger.add_genome('s3', [('s3-1', 'TCGGCCT')])
        merger.add_snp('s1', 's1-1', 8, [('A', 80, 80), ('C', 20, 20)])
        merger.add_snp('s2', 's2-1', 7, [('A', 90, 90), ('C', 10, 10)])
        merger.add_snp('s3', 's3-1', 5, [('C', 70, 70), ('A', 30, 30)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 6)
        self.assertEqual(rows[0].ref, 'A')
        self.assertEqual(rows[0].alt, 'C')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.2')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '0:0.1')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '1:0.7')

    def test_sample_major_allele_not_ref_allele(self):
        # make sure we can invert the allele frequency of the isnv
        # if necessary to match the reference's definition of ref & alt
        merger = VcfMergeRunner([('ref1', 'ATCGCAC')])
        merger.add_genome('s1', [('s1-1', 'ATAGCCC')])
        merger.add_snp('s1', 's1-1', 3, [('C', 10, 10), ('A', 90, 90)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 3)
        self.assertEqual(rows[0].ref, 'C')
        self.assertEqual(rows[0].alt, 'A')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '1:0.9')

    def test_backfill_sample_from_assembly(self):
        # one sample has no isnv, but another does, so we fill it in
        # 100% with its assembly allele, unless it doesn't have one
        # REF C
        # S1  A (isnv)
        # S2  A (consensus, no isnv)
        # S3  N (no consensus, no isnv)
        merger = VcfMergeRunner([('ref1', 'ATCG')])
        merger.add_genome('s1', [('s1-1', 'ATCG')])
        merger.add_genome('s2', [('s2-1', 'ATAG')])
        merger.add_genome('s3', [('s3-1', 'ATNG')])
        merger.add_snp('s1', 's1-1', 3, [('C', 90, 90), ('A', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 3)
        self.assertEqual(rows[0].ref, 'C')
        self.assertEqual(rows[0].alt, 'A')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.1')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '1:1.0')
        self.assertEqual(rows[0][2], '.:.:.:.:.')

    def test_simple_insertions(self):
        # IA, ITCG, etc
        # V-Phaser outputs a position that is just prior to where the new
        # bases come in.  For example: ATCG -> ATiiiCG is at position 2.
        # This is consistent with what the VCF conventional position.
        # V-Phaser outputs an allele that does not include the position.
        # For example: ATCG -> ATAAACG is considered "IAAA" at position 2.
        # This is not the same as the VCF convention, which includes the
        # initial invariant base as part of the allele (it's a T -> TAAA
        # variant at position 2).
        merger = VcfMergeRunner([('ref1', 'ATCG')])
        merger.add_genome('s1', [('s1-1', 'ATCG')])
        merger.add_genome('s2', [('s2-1', 'ATCG')])
        merger.add_genome('s3', [('s3-1', 'ATCG')])
        merger.add_indel('s1', 's1-1', 2, [('', 80, 80), ('AA', 20, 20)])
        merger.add_indel('s2', 's2-1', 2, [('', 90, 90), ('AT', 10, 10)])
        merger.add_indel('s3', 's3-1', 2, [('', 80, 80), ('GCC', 15, 15), ('AA', 5, 5)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 2)
        self.assertEqual(rows[0].ref, 'T')
        self.assertEqual(rows[0].alt, 'TAA,TGCC,TAT')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.2,0.0,0.0')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '0:0.0,0.0,0.1')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '0:0.05,0.15,0.0')

    def test_simple_deletions(self):
        # D1, D2, etc...
        # V-Phaser outputs a position that describes the deleted base and
        # the number of bases deleted (but does not tell you what bases
        # were deleted).  For example: ATCGAT -> ATAT is considered
        # a "D2" deletion at position 3.
        # This is not the same as the VCF convention, which anchors on
        # a preceeding invariant base. The above example is considered
        # to be a TCG -> T variant at position 2.
        merger = VcfMergeRunner([('ref1', 'ATCGAT')])
        merger.add_genome('s1', [('s1-1', 'ATCGAT')])
        merger.add_genome('s2', [('s2-1', 'ATCGAT')])
        merger.add_genome('s3', [('s3-1', 'ATCGAT')])
        merger.add_indel('s1', 's1-1', 3, [('CG', 80, 80), ('', 20, 20)])
        merger.add_indel('s2', 's2-1', 3, [('C', 90, 90), ('', 10, 10)])
        merger.add_indel('s3', 's3-1', 3, [('CG', 80, 80), ('G', 15, 15), ('', 5, 5)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 2)
        self.assertEqual(rows[0].ref, 'TCG')
        self.assertEqual(rows[0].alt, 'TG,T')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '0:0.0,0.2')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '0:0.1,0.0')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '0:0.15,0.05')

    def test_deletion_spans_deletion(self):
        # sample assembly has deletion against reference and isnv deletes even more
        # POS is anchored right before the deletion
        # REF:  ATCGTTCA
        # S1:   ATCG--CA
        # isnv:       x  (position 5, D1)
        merger = VcfMergeRunner([('ref1', 'ATCGTTCA')])
        merger.add_genome('s1', [('s1-1', 'ATCGCA')])
        merger.add_indel('s1', 's1-1', 5, [('C', 90, 90), ('', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'GTTC')
        self.assertEqual(rows[0].alt, 'GC,G')
        self.assertEqual(rows[0][0], '1:0.9,0.1:200:0,1,1:.,1.0,1.0')

    def test_insertion_spans_deletion(self):
        # sample assembly has deletion against reference and isnv inserts back into it
        # POS is anchored right before the deletion
        # REF:  ATCGTTGA
        # S1:   ATCG--GA
        # isnv:     T     (position 4, IT)
        # isnv:     TT    (position 4, ITT)
        # isnv:     TTC   (position 4, ITTC)
        merger = VcfMergeRunner([('ref1', 'ATCGTTCACC')])
        merger.add_genome('s1', [('s1-1', 'ATCGCACC')])
        merger.add_genome('s2', [('s2-1', 'ATCGCACC')])
        merger.add_genome('s3', [('s3-1', 'ATCGCACC')])
        merger.add_indel('s1', 's1-1', 4, [('', 70, 70), ('T', 30, 30)])
        merger.add_indel('s2', 's2-1', 4, [('', 80, 80), ('TT', 20, 20)])
        merger.add_indel('s3', 's3-1', 4, [('', 90, 90), ('TTC', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'GTT')
        self.assertEqual(rows[0].alt, 'G,GT,GTTC')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '1:0.7,0.3,0.0')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '1:0.8,0.0,0.0')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '1:0.9,0.0,0.1')

    def test_snp_within_insertion(self):
        # sample assembly has insertion against reference and isnp modifies it
        # REF:  ATCG--GA
        # S1:   ATCGTTGA
        # isnv:    C
        # isnv:     C
        # isnv:      C
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGTTGACT')])
        merger.add_genome('s2', [('s2-1', 'TCGTTGACT')])
        merger.add_genome('s3', [('s3-1', 'TCGTTGACT')])
        merger.add_snp('s1', 's1-1', 4, [('G', 70, 70), ('C', 30, 30)])
        merger.add_snp('s2', 's2-1', 4, [('T', 80, 80), ('C', 20, 20)])
        merger.add_snp('s3', 's3-1', 5, [('T', 90, 90), ('C', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'G')
        self.assertEqual(rows[0].alt, 'GTT,CTT,GCT,GTC')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '1:0.7,0.3,0.0,0.0')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '1:0.8,0.0,0.2,0.0')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '1:0.9,0.0,0.0,0.1')

    def test_2snps_within_insertion_same_sample(self):
        # Sample assembly has insertion against reference containing two SNPs
        # in the same sample.
        # We don't handle this because we can't see phasing; it should be
        #     handled by the variant caller.
        # REF:  ATCG--GACT
        # S1:   ATCGTTGACT
        # isnv:     C
        # isnv:      C
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGTTGACT')])
        merger.add_snp('s1', 's1-1', 5, [('T', 80, 80), ('C', 20, 20)])
        merger.add_snp('s1', 's1-1', 6, [('T', 90, 90), ('C', 10, 10)])
        self.assertRaises(NotImplementedError, merger.run_and_get_vcf_rows)

    def test_deletion_past_end_of_some_consensus(self):
        # Some sample contains a deletion beyond the end of the consensus
        # sequence of another sample with a SNP. It should skip latter rather
        # than crashing.
        # REF:    ATCGAAC
        # S1:     ATCG--C
        # S1isnv: ATCA--C
        # S2:     ATCT
        # S2isnv: ATCC
        merger = VcfMergeRunner([('ref1', 'ATCGAAC')])
        merger.add_genome('s1', [('s1-1', 'ATCGC')])
        merger.add_snp('s1', 's1-1', 4, [('G', 70, 70), ('A', 30, 30)])
        merger.add_genome('s2', [('s2-1', 'ATCT')])
        merger.add_snp('s2', 's2-1', 4, [('T', 80, 80), ('C', 20, 20)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'GAA')
        # multiple options because the allele frequencies can be the same
        self.assertIn(rows[0].alt, ['C,G,A', 'G,C,A'])
        self.assertIn(rows[0][0], ['2:0.0,0.7,0.3:200:0,0,1,1:.,.,1.0,1.0', '1:0.7,0.0,0.3:200:0,1,0,1:.,1.0,.,1.0'])
        self.assertIn(rows[0][1], ['1:1.0,0.0,0.0:.:.:.', '2:0.0,1.0,0.0:.:.:.'])
        self.assertEqual(rows[1].contig, 'ref1')
        self.assertEqual(rows[1].pos + 1, 7)
        self.assertEqual(rows[1].ref, 'C')
        self.assertEqual(rows[1].alt, 'T')
        self.assertEqual(rows[1][0], '0:0.0:.:.:.')
        self.assertEqual(rows[1][1], '1:0.8:200:1,1:1.0,1.0')

    def test_snp_past_end_of_some_consensus(self):
        # Some sample contains SNP beyond the end of the consensus
        # sequence of another. It should skip latter rather
        # than crashing.
        # REF:    AT
        # S1:     AT
        # S1isnv: AG
        # S2:     A
        merger = VcfMergeRunner([('ref1', 'AT')])
        merger.add_genome('s1', [('s1-1', 'AT')])
        merger.add_snp('s1', 's1-1', 2, [('T', 70, 70), ('G', 30, 30)])
        merger.add_genome('s2', [('s2-1', 'A')])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 2)
        self.assertEqual(rows[0].ref, 'T')
        self.assertEqual(rows[0].alt, 'G')
        self.assertEqual(rows[0][0], '0:0.3:200:1,1:1.0,1.0')
        self.assertEqual(rows[0][1], '.:.:.:.:.')

    def test_deletion_within_insertion(self):
        # sample assembly has insertion against reference and isnv deletes from it
        # REF:  ATCG--GA
        # S1:   ATCGTTGA
        # isnv:     x     (position 5, D1)  s1          => GTG
        # isnv:     xx    (position 5, D2)  s1          => GG
        # isnv:     xxx   (position 5, D3)  s1          => G
        # isnv:      x    (position 6, D1)  s2, s3      => GTG
        # isnv:      xx   (position 6, D2)  s3          => GT
        # isnv:       x   (position 7, D1)  s4          => GTT
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGTTGACT')])
        merger.add_genome('s2', [('s2-1', 'TCGTTGACT')])
        merger.add_genome('s3', [('s3-1', 'TCGTTGACT')])
        merger.add_genome('s4', [('s4-1', 'TCGTTGACT')])
        merger.add_indel('s1', 's1-1', 5, [('TTG', 40, 40), ('TG', 30, 30), ('G', 20, 20), ('', 10, 10)])
        merger.add_indel('s2', 's2-1', 4, [('T', 80, 80), ('', 20, 20)])
        merger.add_indel('s3', 's3-1', 5, [('TG', 85, 85), ('G', 10, 10), ('', 5, 5)])
        merger.add_indel('s4', 's4-1', 6, [('G', 90, 90), ('', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'GG')
        self.assertEqual(rows[0].alt, 'GTTG,GTG,GTT,G,GT')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '1:0.4,0.3,0.0,0.1,0.0')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '1:0.8,0.2,0.0,0.0,0.0')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '1:0.85,0.1,0.0,0.0,0.05')
        self.assertEqual(':'.join(rows[0][3].split(':')[:2]), '1:0.9,0.0,0.1,0.0,0.0')

    def test_insertion_within_insertion(self):
        # sample assembly has insertion against reference and isnv puts even more in
        # REF:  ATCG--GA
        # S1:   ATCGTTGA
        # isnv:    ^      (position 4, IA)
        # isnv:     ^     (position 5, IA)
        # isnv:      ^    (position 6, IA)
        merger = VcfMergeRunner([('ref1', 'ATCGGACT')])
        merger.add_genome('s1', [('s1-1', 'ATCGTTGACT')])
        merger.add_genome('s2', [('s2-1', 'TCGTTGACT')])
        merger.add_genome('s3', [('s3-1', 'TCGTTGACT')])
        merger.add_indel('s1', 's1-1', 4, [('', 70, 70), ('A', 30, 30)])
        merger.add_indel('s2', 's2-1', 4, [('', 80, 80), ('A', 20, 20)])
        merger.add_indel('s3', 's3-1', 5, [('', 90, 90), ('A', 10, 10)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 4)
        self.assertEqual(rows[0].ref, 'G')
        self.assertEqual(rows[0].alt, 'GTT,GATT,GTAT,GTTA')
        self.assertEqual(':'.join(rows[0][0].split(':')[:2]), '1:0.7,0.3,0.0,0.0')
        self.assertEqual(':'.join(rows[0][1].split(':')[:2]), '1:0.8,0.0,0.2,0.0')
        self.assertEqual(':'.join(rows[0][2].split(':')[:2]), '1:0.9,0.0,0.0,0.1')

    def test_indel_collapse(self):
        # vphaser describes insertions and deletions separately
        # test appropriate collapse of coincident insertions and deletions into
        # a single output VCF row
        # isnv:           (position 2, IA)
        # isnv:           (position 3, D1)
        # Note: This is where V-Phaser gets weird. For sites that have both
        # insertions and deletions, the insertions are described on one row
        # and the deletions get described on a separate row. But both rows
        # will have a read count for the majority/consensus allele
        # (called "i" in one and "d" in the other), but the counts for
        # that allele often do not agree between the two rows!
        # So in this scenario, we ought to average them.
        merger = VcfMergeRunner([('ref1', 'ATCG')])
        merger.add_genome('s1', [('s1-1', 'ATCG')])
        merger.add_indel('s1', 's1-1', 2, [('', 40, 40), ('A', 20, 20)])
        merger.add_indel('s1', 's1-1', 3, [('C', 60, 60), ('', 30, 30)])
        rows = merger.run_and_get_vcf_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].contig, 'ref1')
        self.assertEqual(rows[0].pos + 1, 2)
        self.assertEqual(rows[0].ref, 'TC')
        self.assertEqual(rows[0].alt, 'T,TAC')
        self.assertEqual(rows[0][0].split(':')[0], '0')  # s1 is 0.5 TC, 0.3 T, 0.2 TAC
        for actual, expected in zip(rows[0][0].split(':')[1].split(','), [0.3, 0.2]):
            self.assertAlmostEqual(float(actual), expected, places=2)
