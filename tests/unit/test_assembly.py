# Unit tests for assembly.py

__author__ = "dpark@broadinstitute.org"

import assembly
import util.cmd
import util.file
import Bio.SeqIO
import Bio.Data.IUPACData
import unittest
import argparse
import os
import shutil
import tempfile
import argparse
import itertools
import tools.mummer
import tools.novoalign
import tools.picard
from test import TestCaseWithTmp


def makeFasta(seqs, outFasta):
    with open(outFasta, 'wt') as outf:
        for line in util.file.fastaMaker(seqs):
            outf.write(line)


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in assembly.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()



class TestAmbiguityBases(unittest.TestCase):

    def test_non_failure(self):
        ''' Make sure that alleles_to_ambiguity runs without errors for every possible
            combination of inputs.  Check that the output is one-character long and uppercase.
        '''
        bases = ('A', 'C', 'T', 'G')
        for i in range(1, 5):
            for alleles in itertools.permutations(bases, i):
                out = assembly.alleles_to_ambiguity(alleles)
                self.assertEqual(1, len(out))
                self.assertEqual(out, out.upper())


class TestOrderAndOrient(TestCaseWithTmp):
    ''' Test the MUMmer-based order_and_orient command '''

    def test_varicella_big(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.fasta')
        assembly.order_and_orient(
            os.path.join(inDir, 'contigs.hhv3.fasta'),
            os.path.join(inDir, 'ref.hhv3.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_lassa_multisegment(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.lasv.fasta')
        assembly.order_and_orient(
            os.path.join(inDir, 'contigs.lasv.fasta'),
            os.path.join(inDir, 'ref.lasv.fasta'),
            outFasta)
        self.assertEqualContents(outFasta, expected)
        os.unlink(outFasta)
        
    @unittest.skip('promer alignments not implemented for custom scaffolding step')
    def test_lassa_protein(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.lasv.promer.fasta')
        assembly.order_and_orient(
            os.path.join(inDir, 'contigs.lasv.fasta'),
            os.path.join(inDir, 'ref.lasv.fasta'),
            outFasta,
            aligner='promer')
        self.assertEqualContents(outFasta, expected)
        os.unlink(outFasta)
    
    def test_multi_overlap(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.ebov.small.fasta')
        assembly.order_and_orient(
            os.path.join(inDir, 'contigs.ebov.fasta'),
            os.path.join(inDir, 'ref.ebov.small.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))
        

class TestImputeFromReference(TestCaseWithTmp):
    ''' Test the impute_from_reference command (align and modify_contig) '''

    @unittest.skip('requires 10 mins and 16GB RAM')
    def test_varicella_big_muscle(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.muscle.fasta')
        inDirBase = util.file.get_test_input_path()
        assembly.impute_from_reference(
            os.path.join(inDirBase, 'TestOrderAndOrient', 'expected.hhv3.fasta'),
            os.path.join(inDirBase, 'TestOrderAndOrient', 'ref.hhv3.fasta'),
            outFasta,
            minLengthFraction=0.8,
            minUnambig=0.6,
            replaceLength=55,
            newName='HHV3-test')
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_varicella_big_mummer(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.mummer.fasta')
        inDirBase = util.file.get_test_input_path()
        assembly.impute_from_reference(
            os.path.join(inDirBase, 'TestOrderAndOrient', 'expected.hhv3.fasta'),
            os.path.join(inDirBase, 'TestOrderAndOrient', 'ref.hhv3.fasta'),
            outFasta,
            minLengthFraction=0.8,
            minUnambig=0.6,
            replaceLength=55,
            aligner='mummer',
            newName='HHV3-test')
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_small_muscle(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        assembly.impute_from_reference(
            os.path.join(inDir, 'test.pseudo.fasta'),
            os.path.join(inDir, 'ref.sub.ebov.fasta'),
            outFasta,
            minLengthFraction=0.8,
            minUnambig=0.2,
            replaceLength=5,
            newName='test_sub-EBOV.genome',
            aligner='muscle')
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_small_mafft(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        assembly.impute_from_reference(
            os.path.join(inDir, 'test.pseudo.fasta'),
            os.path.join(inDir, 'ref.sub.ebov.fasta'),
            outFasta,
            minLengthFraction=0.8,
            minUnambig=0.2,
            replaceLength=5,
            newName='test_sub-EBOV.genome',
            aligner='mafft')
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_small_mummer(self):
        inDir = util.file.get_test_input_path(self)
        outFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        assembly.impute_from_reference(
            os.path.join(inDir, 'test.pseudo.fasta'),
            os.path.join(inDir, 'ref.sub.ebov.fasta'),
            outFasta,
            minLengthFraction=0.8,
            minUnambig=0.2,
            replaceLength=5,
            newName='test_sub-EBOV.genome',
            aligner='mummer')
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

class TestRefineAssembly(TestCaseWithTmp):
    def test_ebov_refine1(self):
        inDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'impute.ebov.fasta')
        imputeFasta = util.file.mkstempfname('.imputed.fasta')
        refine1Fasta = util.file.mkstempfname('.refine1.fasta')
        shutil.copy(inFasta, imputeFasta)
        tools.picard.CreateSequenceDictionaryTool().execute(imputeFasta)
        tools.novoalign.NovoalignTool().index_fasta(imputeFasta)
        assembly.refine_assembly(
            imputeFasta,
            os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam'),
            refine1Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 30 -x 20 -t 502',
            min_coverage=2,
            threads=4)
        actual = str(Bio.SeqIO.read(refine1Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine1.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)

    def test_ebov_refine2(self):
        inDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'expected.ebov.refine1.fasta')
        refine1Fasta = util.file.mkstempfname('.refine1.fasta')
        refine2Fasta = util.file.mkstempfname('.refine2.fasta')
        shutil.copy(inFasta, refine1Fasta)
        tools.picard.CreateSequenceDictionaryTool().execute(refine1Fasta)
        tools.novoalign.NovoalignTool().index_fasta(refine1Fasta)
        assembly.refine_assembly(
            refine1Fasta,
            os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam'),
            refine2Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 40 -x 20 -t 100',
            min_coverage=3,
            threads=4)
        actual = str(Bio.SeqIO.read(refine2Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine2.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)


class TestMutableSequence(unittest.TestCase):
    ''' Test the MutableSequence class '''

    def test_bad_coords(self):
        self.assertRaises(Exception, assembly.MutableSequence, 'chr', 0, 4)
        self.assertRaises(Exception, assembly.MutableSequence, 'chr', 5, 4)
        self.assertRaises(Exception, assembly.MutableSequence, 'chr', -2, 4)
        self.assertRaises(Exception, assembly.MutableSequence, 'chr', 5, 6, 'G')

    def test_good_coords(self):
        x = assembly.MutableSequence('chr', 1, 5)
        x = assembly.MutableSequence('chr', 5, 5)
        x = assembly.MutableSequence('chr', 100, 2000)
        x = assembly.MutableSequence('chr name with spaces 5 @#$ --', 1, 5)
        x = assembly.MutableSequence('chr', 5, 5, 'A')
        x = assembly.MutableSequence('chr', 5, 6, 'AT')

    def test_modify_one(self):
        x = assembly.MutableSequence('chr', 5, 8, 'ATCG')
        self.assertRaises(Exception, x.modify, 4, 'G')
        self.assertRaises(Exception, x.modify, 9, 'G')
        self.assertEqual(x.emit(), ('chr', 'ATCG'))
        x.modify(5, 'G')
        self.assertEqual(x.emit(), ('chr', 'GTCG'))
        x.modify(6, 'G')
        self.assertEqual(x.emit(), ('chr', 'GGCG'))
        x.modify(7, 'G')
        self.assertEqual(x.emit(), ('chr', 'GGGG'))
        x.modify(8, 'G')
        self.assertEqual(x.emit(), ('chr', 'GGGG'))
        x.modify(6, 'j')
        self.assertEqual(x.emit(), ('chr', 'GjGG'))
        x.modify(8, 'Y')
        self.assertEqual(x.emit(), ('chr', 'GjGY'))

    def test_modify_blank(self):
        x = assembly.MutableSequence('chr', 5, 8)
        self.assertEqual(x.emit(), ('chr', 'NNNN'))
        x.modify(6, 'G')
        self.assertEqual(x.emit(), ('chr', 'NGNN'))

    def test_modify_insertions(self):
        x = assembly.MutableSequence('chr', 5, 8, 'ATCG')
        x.modify(6, 'insert')
        self.assertEqual(x.emit(), ('chr', 'AinsertCG'))
        x.modify(8, 'tail')
        self.assertEqual(x.emit(), ('chr', 'AinsertCtail'))
        x.modify(5, 'headA')
        self.assertEqual(x.emit(), ('chr', 'headAinsertCtail'))

    def test_modify_deletions(self):
        x = assembly.MutableSequence('chr', 5, 8, 'ATCG')
        self.assertRaises(Exception, x.replace, 6, 9, 'AT')
        x.replace(6, 7, 'CT')
        self.assertEqual(x.emit(), ('chr', 'ACTG'))
        x.replace(6, 7, '')
        self.assertEqual(x.emit(), ('chr', 'AG'))
        x.modify(7, 'x')
        self.assertEqual(x.emit(), ('chr', 'AxG'))
        x.modify(6, 'y')
        self.assertEqual(x.emit(), ('chr', 'AyxG'))
        x.replace(7, 7, '123')
        self.assertEqual(x.emit(), ('chr', 'Ay123G'))
        x.modify(7, 'z')
        self.assertEqual(x.emit(), ('chr', 'AyzG'))

    def test_modify_deletions_simple(self):
        x = assembly.MutableSequence('chr', 5, 8, 'ATCG')
        x.replace(6, 7, 'T')
        self.assertEqual(x.emit(), ('chr', 'ATG'))

    def test_modify_deletions_remember(self):
        x = assembly.MutableSequence('chr', 5, 8, 'ATCG')
        x.replace(6, 7, 'T')
        self.assertEqual(x.emit(), ('chr', 'ATG'))
        x.modify(7, 'x')
        self.assertEqual(x.emit(), ('chr', 'ATxG'))
        x.replay_deletions()
        self.assertEqual(x.emit(), ('chr', 'ATG'))


class TestManualSnpCaller(unittest.TestCase):
    ''' Test the vcfrow_parse_and_call_snps method.. lots of edge cases. '''

    def test_missing_dp(self):
        ''' VCF files might contain rows with no calls or any kind of data and that's okay. '''
        row = ['chr10', '105', '.', 'G', '.', '.', '.', '.', 'GT', './.']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])

    def test_dp_inaccurate(self):
        ''' The DP might not equal the sum of the ADs and that's okay apparently. '''
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:5:2,2']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(set(out[0][4]), set(['G', 'A']))
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:2:3,3']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'A']))
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:10:2,0']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [])

    def test_invariant_sites(self):
        ''' Invariant site handling is slightly different in code, so test it specially. '''
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', '0/0:3']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [('LASV.l', 1, 1, 's1', ['T'])])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=0))
        self.assertEqual(out, [('LASV.l', 1, 1, 's1', ['T'])])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', './.']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', './.:10']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [('LASV.l', 1, 1, 's1', ['T'])])

    def test_het_edgecases(self):
        ''' The interplay between min_coverage and major_cutoff is not obvious, here's
            what I understand from Kristian about the desired behavior.
            for min_dp=3:
                3G, 4A, 5C ->  G/A/C
                2G, 3A, 3T -> A/T
                2A, 2T -> no call
                2G, 3C -> C
                2A, 3C, 4T -> T
            for min_dp=2:
                2A, 2T -> A/T
         '''
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:3,4,5,0']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'A', 'C']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,3,0,3']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['A', 'T']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [])
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=2))
        self.assertEqual(set(out[0][4]), set(['A', 'T']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,0,3,0']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out[0][4], ['C'])
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,3,4']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out[0][4], ['T'])

    def test_indels(self):
        ''' Indel handling '''
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,10,1']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['GA']))
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,5,2']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'GA']))
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,5,3']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'GA', 'T']))
        row = ['thecontig', '105000', '.', 'AT', 'A', '.', '.', '.', 'GT:AD', '0/1:2,10']
        out = list(assembly.vcfrow_parse_and_call_snps(row, ['s1']))
        self.assertEqual(out, [('thecontig', 105000, 105001, 's1', ['A'])])

    def test_vcf_to_seqs_indels1(self):
        input = ['thecontig', '5', '.', 'AT', 'A', '.', '.', '.', 'GT:AD', '0/1:2,10']
        actual = assembly.vcf_to_seqs([input], {'thecontig': 10}, ['s1'], min_dp=2)
        actual = list(actual)[0][1].strip('N')
        self.assertEqual(actual, 'A')
        actual = assembly.vcf_to_seqs([input], {'thecontig': 10}, ['s1'], min_dp=2)
        actual = list(actual)[0][1]
        self.assertEqual(actual, 'NNNNANNNN')

    def test_vcf_to_seqs_indels2(self):
        ''' More end-to-end indel handling '''
        myInputDir = util.file.get_test_input_path(self)
        input = os.path.join(myInputDir, 'indel.vcf.gz')
        expected = os.path.join(myInputDir, 'output.fasta')
        chrlens = {'EBOV_2014_G6060.1': 18962}
        samples = ['G6060.1']
        expected = str(Bio.SeqIO.read(expected, 'fasta').seq)
        actual = assembly.vcf_to_seqs(util.file.read_tabfile(input), chrlens, samples, min_dp=2)
        actual = list(actual)[0][1].strip('N')
        self.assertEqual(actual, expected)


class TestDeambigAndTrimFasta(TestCaseWithTmp):
    ''' Test the deambig_fasta and trim_fasta commands. '''

    def run_method(self, inseqs, parser_fun):
        fasta_in = util.file.mkstempfname()
        fasta_out = util.file.mkstempfname()
        makeFasta([(str(i), inseqs[i]) for i in range(len(inseqs))], fasta_in)
        args = parser_fun(argparse.ArgumentParser()).parse_args([fasta_in, fasta_out])
        args.func_main(args)
        return (fasta_in, fasta_out)

    def test_trim_fasta(self):
        ''' Simple test of the trim_fasta command '''
        inseqs = ['NNnnNNnNaslkdfjasdkfNNNN', 'NNNnnN', 'NNN123', 'ATCG']
        expected = ['aslkdfjasdkf', '', '123', 'ATCG']
        expected = dict((str(i), expected[i]) for i in range(len(expected)))
        fasta_in, fasta_out = self.run_method(inseqs, assembly.parser_trim_fasta)
        with open(fasta_out, 'rt') as fa:
            for record in Bio.SeqIO.parse(fa, 'fasta'):
                self.assertIn(record.id, expected)
                self.assertEqual(str(record.seq), expected[record.id])

    def test_deambig_fasta(self):
        ''' Simple test of the deambig_fasta command '''
        table = [(k, v) for k, v in Bio.Data.IUPACData.ambiguous_dna_values.items() if k != 'X']
        keys = [k for k, v in table]
        vals = [set(v) for k, v in table]
        keys = keys + [k.lower() for k in keys]
        vals = vals + vals
        inseq = ''.join(keys)
        fasta_in, fasta_out = self.run_method([inseq], assembly.parser_deambig_fasta)
        with open(fasta_out, 'rt') as fa:
            for rec in Bio.SeqIO.parse(fa, 'fasta'):
                self.assertEqual(rec.id, '0')
                outseq = str(rec.seq)
                for i in range(len(outseq)):
                    self.assertIn(outseq[i], vals[i])


class TestContigChooser(unittest.TestCase):
    ''' Test the contig_chooser heuristic used by our MUMmer-based custom scaffolder. '''

    def test_no_seqs(self):
        for test_len in (7,2,228,52):
            actual = tools.mummer.contig_chooser([], test_len)
            self.assertEqual(actual, 'N' * test_len)

    def test_one_seq(self):
        for test_seq in ('A', '', 'GACTGATG', 'non-biological :characters!'):
            actual = tools.mummer.contig_chooser([test_seq], 90)
            self.assertEqual(actual, test_seq)
    
    def test_most_popular_seq(self):
        alt_seqs = ['AA', 'aa', 'GGA', 'T', 'GGA']
        expected = 'GGA'
        actual = tools.mummer.contig_chooser(alt_seqs, 2)
        self.assertEqual(actual, expected)

    def test_most_popular_seq_len(self):
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC', 'aa']
        actual = tools.mummer.contig_chooser(alt_seqs, 2)
        self.assertEqual(actual, 'aa')
        actual = tools.mummer.contig_chooser(alt_seqs, 3)
        self.assertEqual(actual, 'GGA')
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC']
        actual = tools.mummer.contig_chooser(alt_seqs, 20)
        self.assertEqual(actual, 'GGA')
        actual = tools.mummer.contig_chooser(alt_seqs, 1)
        self.assertEqual(actual, 'GGA')
    
    def test_same_as_ref_len(self):
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC', 'aa']
        actual = tools.mummer.contig_chooser(alt_seqs, 1)
        self.assertEqual(actual, 'T')


