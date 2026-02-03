# Unit tests for viral_ngs.assembly.py

__author__ = "dpark@broadinstitute.org"

import viral_ngs.assembly
import viral_ngs.core.cmd
import viral_ngs.core.file
import Bio.SeqIO
import Bio.Data.IUPACData
import unittest
import argparse
import csv
import glob
import os
import os.path
import shutil
import subprocess
import tempfile
import argparse
import itertools
import pytest
import viral_ngs.assemble.mummer
import viral_ngs.assemble.skani
import viral_ngs.core.minimap2
import viral_ngs.core.novoalign
import viral_ngs.core.picard
import viral_ngs.core.samtools
from tests import TestCaseWithTmp, _CPUS


def makeFasta(seqs, outFasta):
    with open(outFasta, 'wt') as outf:
        for line in viral_ngs.core.file.fastaMaker(seqs):
            outf.write(line)


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in viral_ngs.assembly.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestRefineAssemble(TestCaseWithTmp):
    ''' Test edge cases of the de novo assembly pipeline '''

    def test_empty_input_bam_assembly(self):
        novoalign = viral_ngs.core.novoalign.NovoalignTool()
        novoalign.install()

        # prep inputs
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta')
        inFasta = viral_ngs.core.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, inFasta)
        novoalign.index_fasta(inFasta)

        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')

        outFasta = viral_ngs.core.file.mkstempfname('.refined.fasta')

        # run refine_assembly
        args = [inFasta, inBam, outFasta, "--chr_names", 'G5012.3', "--min_coverage", '3', "--novo_params",
                "-r Random -l 30 -g 40 -x 20 -t 502 -c {}".format(_CPUS)]
        args = viral_ngs.assembly.parser_refine_assembly(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # the expected output is an empty fasta file
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) == 0)

    def test_aligned_empty_input_bam_assembly(self):
        mm2 = viral_ngs.core.minimap2.Minimap2()
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        with viral_ngs.core.file.tempfname('.ref.fasta') as inFasta:
            with viral_ngs.core.file.tempfname('.bam') as inBam:
                with viral_ngs.core.file.tempfname('.refined.fasta') as outFasta:
                    with viral_ngs.core.file.tempfname('.vcf.gz') as outVcf:
                        shutil.copyfile(os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta'), inFasta)
                        mm2.align_bam(os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam'), inFasta, inBam,
                            options=['-x', 'sr'])
                        samtools.index(inBam)

                        # run refine_assembly
                        args = [inFasta, inBam, outFasta, "--already_realigned_bam", inBam, '--outVcf', outVcf]
                        args = viral_ngs.assembly.parser_refine_assembly(argparse.ArgumentParser()).parse_args(args)
                        args.func_main(args)

                        # the expected output is an empty fasta file
                        self.assertTrue(os.path.isfile(outFasta))
                        self.assertTrue(os.path.getsize(outFasta) == 0)
                        self.assertTrue(os.path.isfile(outVcf))
                        self.assertTrue(os.path.getsize(outVcf) > 0)

    def test_empty_input_fasta_assembly(self):
        novoalign = viral_ngs.core.novoalign.NovoalignTool()
        novoalign.install()

        # make the input fasta empty
        inFasta = viral_ngs.core.file.mkstempfname('.input.fasta')
        viral_ngs.core.file.touch(inFasta)
        novoalign.index_fasta(inFasta)

        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')

        outFasta = viral_ngs.core.file.mkstempfname('.refined.fasta')

        # run refine_assembly
        args = [inFasta, inBam, outFasta, "--chr_names", 'G5012.3', "--min_coverage", '3', "--novo_params",
                "-r Random -l 30 -g 40 -x 20 -t 502 -c {}".format(_CPUS)]
        args = viral_ngs.assembly.parser_refine_assembly(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # the expected output is an empty fasta file
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) == 0)


    def test_empty_input_succeed(self):
        novoalign = viral_ngs.core.novoalign.NovoalignTool()
        novoalign.install()

        # make the input fasta empty
        inFasta = viral_ngs.core.file.mkstempfname('.input.fasta')
        viral_ngs.core.file.touch(inFasta)
        novoalign.index_fasta(inFasta)

        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')

        outFasta = viral_ngs.core.file.mkstempfname('.refined.fasta')

        # run refine_assembly
        args = [inFasta, inBam, outFasta, "--chr_names", 'G5012.3', "--min_coverage", '3', "--novo_params",
                "-r Random -l 30 -g 40 -x 20 -t 502 -c {}".format(_CPUS)]
        args = viral_ngs.assembly.parser_refine_assembly(argparse.ArgumentParser()).parse_args(args)
        print(args)
        args.func_main(args)

        # the expected output is an empty fasta file
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) == 0)


class TestAssembleSpades(TestCaseWithTmp):
    ''' Test the assemble_spades command (no validation of output) '''

    def test_assembly(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inBam = os.path.join(inDir, '..', 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'clipDb.fasta')
        with viral_ngs.core.file.tempfname('.fasta') as outFasta:
            viral_ngs.assembly.assemble_spades(in_bam=inBam, clip_db=clipDb, min_contig_len=180, out_fasta=outFasta)
            self.assertGreater(os.path.getsize(outFasta), 0)
            contig_lens = list(sorted(len(seq.seq) for seq in Bio.SeqIO.parse(outFasta, 'fasta')))
            #import sys
            #print('test_assembly_contigs_lens:', contig_lens, file=sys.stderr)
            self.assertEqual(contig_lens, [200, 201, 210, 222, 243, 244, 268, 294, 300, 328, 348, 376, 381, 430])

    def test_assembly_with_previously_assembled_contigs(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inBam = os.path.join(inDir, '..', 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'clipDb.fasta')
        previously_assembled_contigs = os.path.join(inDir, 'trinity_contigs.fasta')
        with viral_ngs.core.file.tempfname('.fasta') as outFasta:
            viral_ngs.assembly.assemble_spades(in_bam=inBam, clip_db=clipDb, contigs_untrusted=previously_assembled_contigs,
                                     out_fasta=outFasta, mem_limit_gb=1)
            self.assertGreater(os.path.getsize(outFasta), 0)
            contig_lens = list(sorted(len(seq.seq) for seq in Bio.SeqIO.parse(outFasta, 'fasta')))
            #import sys
            #print('test_assembly_with_previously_assembled_contigs_contigs_lens:', contig_lens, file=sys.stderr)
            self.assertEqual(contig_lens, [200, 201, 210, 222, 243, 244, 268, 294, 300, 328, 348, 376, 381, 430])

    def test_empty_input_succeed(self):
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'empty.bam')
        clipDb = os.path.join(inDir, 'clipDb.fasta')
        with viral_ngs.core.file.tempfname('fasta') as outFasta:
            viral_ngs.assembly.assemble_spades(in_bam=inBam, clip_db=clipDb, out_fasta=outFasta)
            self.assertEqual(os.path.getsize(outFasta), 0)

    def test_always_succeed(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inBam = os.path.join(inDir, '..', 'almost-empty.bam')
        clipDb = os.path.join(inDir, 'clipDb.fasta')
        with viral_ngs.core.file.tempfname('.fasta') as outFasta:
            viral_ngs.assembly.assemble_spades(in_bam=inBam, clip_db=clipDb, out_fasta=outFasta, spades_opts='--bad-option',
                                     always_succeed=True)
            self.assertEqual(os.path.getsize(outFasta), 0)


class TestAmbiguityBases(unittest.TestCase):

    def test_non_failure(self):
        ''' Make sure that alleles_to_ambiguity runs without errors for every possible
            combination of inputs.  Check that the output is one-character long and uppercase.
        '''
        bases = ('A', 'C', 'T', 'G')
        for i in range(1, 5):
            for alleles in itertools.permutations(bases, i):
                out = viral_ngs.assembly.alleles_to_ambiguity(alleles)
                self.assertEqual(1, len(out))
                self.assertEqual(out, out.upper())


class TestUndirectedGraph(unittest.TestCase):
    def test_simple(self):
        g = viral_ngs.assemble.skani.UndirectedGraph()
        g.add_edge('a', 'b')
        g.add_edge('a', 'c')
        g.add_edge('b', 'd')
        actual = list(sorted(g.get_clusters()))
        self.assertEqual(actual, [{'a', 'b', 'c', 'd'}])

    def test_disconnected(self):
        g = viral_ngs.assemble.skani.UndirectedGraph()
        g.add_edge('a', 'b')
        g.add_edge('c', 'd')
        actual = list(sorted(g.get_clusters()))
        self.assertEqual(actual, [{'a', 'b'}, {'c', 'd'}])

    def test_both(self):
        g = viral_ngs.assemble.skani.UndirectedGraph()
        g.add_edge(1, 2)
        g.add_edge(11,12)
        g.add_edge(18,15)
        g.add_node(12)
        g.add_node(22)
        g.add_node(55)
        g.add_edge(25,22)
        g.add_edge(7,2)
        g.add_edge(12,18)
        actual = list(sorted(g.get_clusters()))
        self.assertEqual(actual, [{1, 2, 7}, {11, 12, 15, 18}, {22, 25}, {55}])


class TestOrderAndOrient(TestCaseWithTmp):
    ''' Test the MUMmer-based order_and_orient command '''

    def test_varicella_big(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.hhv3.fasta'),
            os.path.join(inDir, 'ref.hhv3.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_lassa_multisegment(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.lasv.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.lasv.fasta'),
            os.path.join(viral_ngs.core.file.get_test_input_path(), 'ref.lasv.fasta'),
            outFasta)
        self.assertEqualContents(outFasta, expected)
        os.unlink(outFasta)

    def test_lassa_multisegment_refsel(self):
        with viral_ngs.core.file.tempfnames(('.out.fasta', '.out_ref.fasta', '.stats.tsv')) \
             as (outFasta, outReference, outStats):
            contigs, expected, expectedStats = self.inputs('contigs.lasv.fasta', 
                                                           'expected.lasv.fasta', 
                                                           'expected.refsel.lasv.stats.tsv')
            refs = [self.input('ref.lasv.{}.fasta'.format(strain))
                    for strain in ('josiah', 'pinneo', 'KGH_G502', 'BNI_Nig08_A19', 'nomatch')]
            viral_ngs.assembly.order_and_orient(contigs, refs, outFasta,
                                      outReference=outReference, outStats=outStats)
            self.assertEqualContents(outFasta, expected)
            self.assertEqualFasta(outReference, refs[0])
            self.assertEqualContents(outStats, expectedStats)

    def test_influenza_multisegment(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.influenza.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.influenza.fasta'),
            os.path.join(inDir, 'ref.influenza.fasta'),
            outFasta)
        self.assertEqualContents(outFasta, expected)
        os.unlink(outFasta)

    def test_ebov_palindrome(self):
        # this tests a scenario where show-aligns has more alignments than show-tiling
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.ebov.doublehit.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.ebov.doublehit.fasta'),
            os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_ebov_palindrome_refsel(self):
        # this tests a scenario where show-aligns has more alignments than show-tiling
        with viral_ngs.core.file.tempfnames(('.out.fasta', '.stats.tsv')) as (outFasta, outStats):
            contigs, expected, expectedStats = self.inputs('contigs.ebov.doublehit.fasta',
                                                                 'expected.ebov.doublehit.fasta',
                                                                 'expected.refsel.ebov.stats.tsv')
            refs = self.inputs('ref.ebov.lbr.fasta','ref.ebov.sle.fasta','ref.ebov.gin.fasta')
            viral_ngs.assembly.order_and_orient(contigs, refs, outFasta, outStats=outStats)
            self.assertEqualFastaSeqs(outFasta, expected)
            self.assertEqualContents(outStats, expectedStats)

    def test_hiv_wraparound(self):
        # this tests a misassembly from Trinity and checks that we still use some of the contig
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hiv.wrapped.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.hiv.wrapped.fasta'),
            os.path.join(inDir, 'ref.hiv.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_alternate_contigs(self):
        # this tests that --outAlternateContigs works as expected
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        altFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hiv.big_indel.fasta')
        expectedAlt = os.path.join(inDir, 'expected.hiv.big_indel.alternates.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.hiv.big_indel.fasta'),
            os.path.join(inDir, 'ref.hiv.fasta'),
            outFasta,
            outAlternateContigs=altFasta)
        self.assertEqualContents(outFasta, expected)
        self.assertEqualContents(altFasta, expectedAlt)

    @unittest.skip('promer alignments not implemented for custom scaffolding step')
    def test_lassa_protein(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.lasv.promer.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.lasv.fasta'),
            os.path.join(viral_ngs.core.file.get_test_input_path(), 'ref.lasv.fasta'),
            outFasta,
            aligner='promer')
        self.assertEqualContents(outFasta, expected)
        os.unlink(outFasta)

    def test_multi_overlap(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.ebov.small.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.ebov.fasta'),
            os.path.join(inDir, 'ref.ebov.small.fasta'),
            outFasta)
        self.assertEqual(
            str(Bio.SeqIO.read(outFasta, 'fasta').seq),
            str(Bio.SeqIO.read(expected, 'fasta').seq))

    def test_ambig_align(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        contigs_gz = os.path.join(inDir, 'contigs.lasv.ambig.fasta.gz')
        contigs = viral_ngs.core.file.mkstempfname('.fasta')
        with viral_ngs.core.file.open_or_gzopen(contigs_gz, 'rb') as f_in:
            with open(contigs, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        expected = os.path.join(inDir, 'expected.lasv.ambig.fasta')
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.order_and_orient(
            contigs,
            os.path.join(inDir, 'ref.lasv.ISTH2376.fasta'),
            outFasta)
        def get_seqs(fasta):
            return [str(s.seq) for s in Bio.SeqIO.parse(fasta, 'fasta')]
        self.assertEqual(get_seqs(outFasta), get_seqs(expected))

    def test_ambig_align_ebov(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        contigs_gz = os.path.join(inDir, 'contigs.ebov.ambig.fasta.gz')
        contigs = viral_ngs.core.file.mkstempfname('.fasta')
        with viral_ngs.core.file.open_or_gzopen(contigs_gz, 'rb') as f_in:
            with open(contigs, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        expected = os.path.join(inDir, 'expected.ebov.ambig.fasta')
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.order_and_orient(
            contigs,
            os.path.join(inDir, 'ref.ebov.makona_C15.fasta'),
            outFasta)
        def get_seqs(fasta):
            return [str(s.seq) for s in Bio.SeqIO.parse(fasta, 'fasta')]
        self.assertEqual(get_seqs(outFasta), get_seqs(expected))

    def test_obscure_mummer3_bug(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.lasv.bug.fasta')
        # under mummer3, this fails with a weird error in nucmer,
        # causing a CalledProcessError. We want this to succeed
        # nucmer, but fail later on due to IncompleteAssemblyError
        self.assertRaises(viral_ngs.assembly.IncompleteAssemblyError,
            viral_ngs.assembly.order_and_orient,
            os.path.join(inDir, 'contig.mummer3_fail_lasv.fasta'),
            os.path.join(inDir, 'ref.lasv.ISTH2376.fasta'),
            outFasta)

    def test_not_all_segments_fail(self):
        # IncompleteAssemblyError is thrown when only some but not all segments are recovered
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        self.assertRaises(viral_ngs.assembly.IncompleteAssemblyError,
            viral_ngs.assembly.order_and_orient,
            os.path.join(inDir, 'contigs.lasv.one_small.fasta'),
            os.path.join(inDir, 'ref.lasv.ISTH2376.fasta'),
            outFasta)

    def test_not_all_segments_succeed(self):
        # ... unless we specifically allow for partial outputs
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.order_and_orient(
            os.path.join(inDir, 'contigs.lasv.one_small.fasta'),
            os.path.join(inDir, 'ref.lasv.ISTH2376.fasta'),
            outFasta,
            allow_incomplete_output=True)

    def test_empty_output_succeed(self):
        # a completely empty output should be possible if we allow it
        inDir = viral_ngs.core.file.get_test_input_path(self)
        emptyFasta = os.path.join(inDir, '..', 'empty.fasta')
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.order_and_orient(
            emptyFasta,
            os.path.join(inDir, 'ref.influenza.fasta'),
            outFasta,
            allow_incomplete_output=True)
        self.assertEqualContents(outFasta, emptyFasta)

class TestGap2Seq(TestCaseWithTmp):
    '''Test gap-filling tool Gap2Seq'''

    def test_gapfill(self):
        join = os.path.join
        inDir = viral_ngs.core.file.get_test_input_path()
        in_scaffold = join(inDir, 'TestOrderAndOrient', 'expected.ebov.doublehit.fasta')
        with viral_ngs.core.file.tempfname(suffix='.filled.fasta') as filled:
            viral_ngs.assembly.gapfill_gap2seq(in_scaffold=in_scaffold,
                                     in_bam=join(inDir, 'G5012.3.testreads.bam'),
                                     out_scaffold=filled, random_seed=23923937, threads=1)
            shutil.copyfile(filled, '/tmp/filled.fasta')
            self.assertEqualContents(filled, join(inDir, 'TestGap2Seq', 'expected.ebov.doublehit.gapfill.fasta'))

    def test_empty_fasta_input(self):
        inDir = viral_ngs.core.file.get_test_input_path()
        empty_fasta = os.path.join(inDir, 'empty.fasta')
        out_fasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.gapfill_gap2seq(in_scaffold=empty_fasta,
                                    in_bam=os.path.join(inDir, 'G5012.3.testreads.bam'),
                                    out_scaffold=out_fasta, random_seed=23923937, threads=1)
        self.assertEqualContents(out_fasta, empty_fasta)


class TestImputeFromReference(TestCaseWithTmp):
    ''' Test the impute_from_reference command (align and modify_contig) '''

    @unittest.skip('requires 10 mins and 16GB RAM')
    def test_varicella_big_muscle(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.muscle.fasta')
        inDirBase = viral_ngs.core.file.get_test_input_path()
        viral_ngs.assembly.impute_from_reference(
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
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.hhv3.mummer.fasta')
        inDirBase = viral_ngs.core.file.get_test_input_path()
        viral_ngs.assembly.impute_from_reference(
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
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        viral_ngs.assembly.impute_from_reference(
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
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        viral_ngs.assembly.impute_from_reference(
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
        inDir = viral_ngs.core.file.get_test_input_path(self)
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(inDir, 'expected.sub.ebov.impute.fasta')
        viral_ngs.assembly.impute_from_reference(
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

    def test_empty_fasta_input(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        empty_fasta = os.path.join(inDir, '..', 'empty.fasta')
        outFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.impute_from_reference(
            empty_fasta,
            os.path.join(inDir, 'ref.sub.ebov.fasta'),
            outFasta,
            minLengthFraction=0.0,
            minUnambig=0.0,
            replaceLength=5,
            newName='test_sub-EBOV.genome',
            aligner='mummer')
        self.assertEqualContents(outFasta, empty_fasta)

class TestSkaniReferenceSelection(TestCaseWithTmp):
    ''' Test Skani-based reference selection '''

    def test_skani_contigs_to_refs(self):
        '''
            Test the skani_contigs_to_refs function.
            Test inputs include LASV MAGs/contigs against various EBOV and LASV references.
            The only references that should hit are the LASV Josiah and KGH_G502 references.
            Additionally, skani should identify them as being from the same cluster.
            No EBOV references should be selected.
        '''

        inDir = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestOrderAndOrient')
        with viral_ngs.core.file.tempfnames(('.skani.dist.out', '.skani.dist.filtered', '.clusters.filtered')) \
            as (out_skani_dist, out_skani_dist_filtered, out_clusters_filtered):
            contigs = os.path.join(inDir, 'contigs.lasv.fasta')
            refs =  [os.path.join(inDir, 'ref.lasv.{}.fasta'.format(strain))
                     for strain in ('josiah', 'pinneo', 'KGH_G502', 'BNI_Nig08_A19', 'nomatch')] + \
                    [os.path.join(inDir, 'ref.ebov.{}.fasta'.format(strain))
                      for strain in ('lbr', 'sle', 'gin')]

            viral_ngs.assembly.skani_contigs_to_refs(contigs, refs, out_skani_dist, out_skani_dist_filtered, out_clusters_filtered, threads=1)

            with open(out_clusters_filtered, 'r') as inf:
                clusters = inf.readlines()
            self.assertEqual(len(clusters), 1)
            actual_cluster = set([os.path.basename(f) for f in clusters[0].strip().split('\t')])
            expected_cluster = set(['ref.lasv.{}.fasta'.format(strain) for strain in ('josiah', 'pinneo', 'KGH_G502', 'BNI_Nig08_A19')])
            self.assertEqual(actual_cluster, expected_cluster)

    def test_skani_no_big_contigs(self):
        '''
            Test that skani_dist tolerates empty (or practically empty) input query fasta
        '''

        inDir = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestOrderAndOrient')
        with viral_ngs.core.file.tempfnames(('.skani.dist.out', '.skani.dist.filtered', '.clusters.filtered')) \
            as (out_skani_dist, out_skani_dist_filtered, out_clusters_filtered):
            contigs = os.path.join(inDir, 'contigs.lasv.one_small.fasta')
            refs =  [os.path.join(inDir, 'ref.ebov.{}.fasta'.format(strain))
                      for strain in ('lbr', 'sle', 'gin')]

            viral_ngs.assembly.skani_contigs_to_refs(contigs, refs, out_skani_dist, out_skani_dist_filtered, out_clusters_filtered, threads=1)

            with open(out_clusters_filtered, 'r') as inf:
                clusters = inf.readlines()
            self.assertEqual(len(clusters), 0)
            with open(out_skani_dist, 'r') as inf:
                lines = inf.readlines()
            self.assertEqual(len(lines), 1)

    def test_skani_no_matches(self):
        '''
            Test that skani_dist tolerates empty outputs (no matches)
        '''

        inDir = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestOrderAndOrient')
        with viral_ngs.core.file.tempfnames(('.skani.dist.out', '.skani.dist.filtered', '.clusters.filtered')) \
            as (out_skani_dist, out_skani_dist_filtered, out_clusters_filtered):
            contigs = os.path.join(inDir, 'contigs.lasv.fasta')
            refs =  [os.path.join(inDir, 'ref.ebov.{}.fasta'.format(strain))
                      for strain in ('lbr', 'sle', 'gin')]

            viral_ngs.assembly.skani_contigs_to_refs(contigs, refs, out_skani_dist, out_skani_dist_filtered, out_clusters_filtered, threads=1)

            with open(out_clusters_filtered, 'r') as inf:
                clusters = inf.readlines()
            self.assertEqual(len(clusters), 0)
            with open(out_skani_dist, 'r') as inf:
                lines = inf.readlines()
            self.assertEqual(len(lines), 1)

    def test_output_sorted_by_product(self):
        '''
            Test that skani.find_closest_reference output tsv is sorted by the product
            of ANI * Total_bases_covered (and not just the default ANI-based sort order)
        '''
        inDir = viral_ngs.core.file.get_test_input_path(self)
        with viral_ngs.core.file.tempfnames(('.skani.dist.out', '.skani.dist.filtered', '.clusters.filtered')) \
            as (out_skani_dist, out_skani_dist_filtered, out_clusters_filtered):
            contigs = os.path.join(inDir, 'USA-MA-Broad_BWH-19947-2023.l000013249603_C5.HTKJ7DRX3.1.acellular.dedup.assembly1-spades.fasta')
            refs =  glob.glob(os.path.join(inDir, 'RVA*.fa'))

            viral_ngs.assembly.skani_contigs_to_refs(contigs, refs, out_skani_dist, out_skani_dist_filtered, out_clusters_filtered, threads=1)

            with open(out_clusters_filtered, 'r') as inf:
                clusters = inf.readlines()
            self.assertEqual(len(clusters), 1)
            with open(out_skani_dist, 'r') as inf:
                hits = list(row['Ref_name'].split()[0].split('_')[1] for row in csv.DictReader(inf, delimiter='\t'))
            self.assertEqual(len(hits), 10)
            expected = ['FJ445177.1','L24917.1','FJ445116.1','DQ473501.1','DQ473496.1','DQ473507.1','DQ473499.1','FJ445140.1','DQ473498.1','GQ223229.1']
            self.assertEqual(hits, expected)

    def test_sort_skani_table_empty_file(self):
        """Test _sort_skani_table_by_product handles completely empty input file"""
        skani = viral_ngs.assemble.skani.SkaniTool()
        with viral_ngs.core.file.tempfnames(('.empty.tsv', '.out.tsv')) as (in_tsv, out_tsv):
            # Create empty file (0 bytes)
            open(in_tsv, 'w').close()

            # Should not raise exception
            skani._sort_skani_table_by_product(in_tsv, out_tsv)

            # Output should exist and have header only (not be empty)
            self.assertTrue(os.path.exists(out_tsv))
            with open(out_tsv) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith('Ref_file'))

    def test_sort_skani_table_header_only(self):
        """Test _sort_skani_table_by_product handles header-only input file"""
        skani = viral_ngs.assemble.skani.SkaniTool()
        header = 'Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\t' \
                 'Ref_name\tQuery_name\tNum_ref_contigs\tNum_query_contigs\tANI_5_percentile\t' \
                 'ANI_95_percentile\tStandard_deviation\tRef_90_ctg_len\tRef_50_ctg_len\t' \
                 'Ref_10_ctg_len\tQuery_90_ctg_len\tQuery_50_ctg_len\tQuery_10_ctg_len\t' \
                 'Avg_chain_len\tTotal_bases_covered\n'

        with viral_ngs.core.file.tempfnames(('.header.tsv', '.out.tsv')) as (in_tsv, out_tsv):
            with open(in_tsv, 'w') as f:
                f.write(header)

            skani._sort_skani_table_by_product(in_tsv, out_tsv)

            # Output should have header only
            with open(out_tsv) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith('Ref_file'))


class TestMutableSequence(unittest.TestCase):
    ''' Test the MutableSequence class '''

    def test_bad_coords(self):
        self.assertRaises(Exception, viral_ngs.assembly.MutableSequence, 'chr', 0, 4)
        self.assertRaises(Exception, viral_ngs.assembly.MutableSequence, 'chr', 5, 4)
        self.assertRaises(Exception, viral_ngs.assembly.MutableSequence, 'chr', -2, 4)
        self.assertRaises(Exception, viral_ngs.assembly.MutableSequence, 'chr', 5, 6, 'G')

    def test_good_coords(self):
        x = viral_ngs.assembly.MutableSequence('chr', 1, 5)
        x = viral_ngs.assembly.MutableSequence('chr', 5, 5)
        x = viral_ngs.assembly.MutableSequence('chr', 100, 2000)
        x = viral_ngs.assembly.MutableSequence('chr name with spaces 5 @#$ --', 1, 5)
        x = viral_ngs.assembly.MutableSequence('chr', 5, 5, 'A')
        x = viral_ngs.assembly.MutableSequence('chr', 5, 6, 'AT')

    def test_modify_one(self):
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8, 'ATCG')
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
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8)
        self.assertEqual(x.emit(), ('chr', 'NNNN'))
        x.modify(6, 'G')
        self.assertEqual(x.emit(), ('chr', 'NGNN'))

    def test_modify_insertions(self):
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8, 'ATCG')
        x.modify(6, 'insert')
        self.assertEqual(x.emit(), ('chr', 'AinsertCG'))
        x.modify(8, 'tail')
        self.assertEqual(x.emit(), ('chr', 'AinsertCtail'))
        x.modify(5, 'headA')
        self.assertEqual(x.emit(), ('chr', 'headAinsertCtail'))

    def test_modify_deletions(self):
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8, 'ATCG')
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
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8, 'ATCG')
        x.replace(6, 7, 'T')
        self.assertEqual(x.emit(), ('chr', 'ATG'))

    def test_modify_deletions_remember(self):
        x = viral_ngs.assembly.MutableSequence('chr', 5, 8, 'ATCG')
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
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])

    def test_dp_inaccurate(self):
        ''' The DP might not equal the sum of the ADs and that's okay apparently. '''
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:5:2,2']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(set(out[0][4]), set(['G', 'A']))
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:2:3,3']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'A']))
        row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:10:2,0']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [])

    def test_invariant_sites(self):
        ''' Invariant site handling is slightly different in code, so test it specially. '''
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', '0/0:3']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [('LASV.l', 1, 1, 's1', ['T'])])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=0))
        self.assertEqual(out, [('LASV.l', 1, 1, 's1', ['T'])])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', './.']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
        self.assertEqual(out, [])
        row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', './.:10']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
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
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'A', 'C']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,3,0,3']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['A', 'T']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out, [])
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=2))
        self.assertEqual(set(out[0][4]), set(['A', 'T']))
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,0,3,0']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out[0][4], ['C'])
        row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,3,4']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(out[0][4], ['T'])

    def test_indels(self):
        ''' Indel handling '''
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,10,1']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['GA']))
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,5,2']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'GA']))
        row = ['thecontig', '105000', '.', 'G', 'GA,T', '.', '.', '.', 'GT:AD', '0/1:5,5,3']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
        self.assertEqual(set(out[0][4]), set(['G', 'GA', 'T']))
        row = ['thecontig', '105000', '.', 'AT', 'A', '.', '.', '.', 'GT:AD', '0/1:2,10']
        out = list(viral_ngs.assembly.vcfrow_parse_and_call_snps(row, ['s1']))
        self.assertEqual(out, [('thecontig', 105000, 105001, 's1', ['A'])])

    def test_vcf_to_seqs_indels1(self):
        input = ['thecontig', '5', '.', 'AT', 'A', '.', '.', '.', 'GT:AD', '0/1:2,10']
        actual = viral_ngs.assembly.vcf_to_seqs([input], {'thecontig': 10}, ['s1'], min_dp=2)
        actual = list(actual)[0][1].strip('N')
        self.assertEqual(actual, 'A')
        actual = viral_ngs.assembly.vcf_to_seqs([input], {'thecontig': 10}, ['s1'], min_dp=2)
        actual = list(actual)[0][1]
        self.assertEqual(actual, 'NNNNANNNN')

    def test_vcf_to_seqs_indels2(self):
        ''' More end-to-end indel handling '''
        myInputDir = viral_ngs.core.file.get_test_input_path(self)
        input = os.path.join(myInputDir, 'indel.vcf.gz')
        expected = os.path.join(myInputDir, 'output.fasta')
        chrlens = {'EBOV_2014_G6060.1': 18962}
        samples = ['G6060.1']
        expected = str(Bio.SeqIO.read(expected, 'fasta').seq)
        actual = viral_ngs.assembly.vcf_to_seqs(viral_ngs.core.file.read_tabfile(input), chrlens, samples, min_dp=2)
        actual = list(actual)[0][1].strip('N')
        self.assertEqual(actual, expected)


class TestDeambigAndTrimFasta(TestCaseWithTmp):
    ''' Test the deambig_fasta and trim_fasta commands. '''

    def run_method(self, inseqs, parser_fun):
        fasta_in = viral_ngs.core.file.mkstempfname()
        fasta_out = viral_ngs.core.file.mkstempfname()
        makeFasta([(str(i), inseqs[i]) for i in range(len(inseqs))], fasta_in)
        args = parser_fun(argparse.ArgumentParser()).parse_args([fasta_in, fasta_out])
        args.func_main(args)
        return (fasta_in, fasta_out)

    def test_trim_fasta(self):
        ''' Simple test of the trim_fasta command '''
        inseqs = ['NNnnNNnNaslkdfjasdkfNNNN', 'NNNnnN', 'NNN123', 'ATCG']
        expected = ['aslkdfjasdkf', '', '123', 'ATCG']
        expected = dict((str(i), expected[i]) for i in range(len(expected)))
        fasta_in, fasta_out = self.run_method(inseqs, viral_ngs.assembly.parser_trim_fasta)
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
        fasta_in, fasta_out = self.run_method([inseq], viral_ngs.assembly.parser_deambig_fasta)
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
            actual = viral_ngs.assemble.mummer.contig_chooser([], test_len)
            self.assertEqual(actual, ['N' * test_len])

    def test_one_seq(self):
        for test_seq in ('A', '', 'GACTGATG', 'non-biological :characters!'):
            actual = viral_ngs.assemble.mummer.contig_chooser([test_seq], 90)
            self.assertEqual(actual, [test_seq])

    def test_most_popular_seq(self):
        alt_seqs = ['AA', 'aa', 'GGA', 'T', 'GGA']
        expected_choice = 'GGA'
        expected_alts = set(('AA', 'aa', 'T'))
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 2)
        self.assertEqual(actual[0], expected_choice)
        self.assertEqual(set(actual[1:]), expected_alts)

    def test_most_popular_seq_len(self):
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC', 'aa']
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 2)
        self.assertEqual(actual[0], 'aa')
        self.assertEqual(set(actual[1:]), set(('AA', 'GGA', 'T', 'GGC')))
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 3)
        self.assertEqual(actual[0], 'GGA')
        self.assertEqual(set(actual[1:]), set(('AA', 'aa', 'T', 'GGC')))
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC']
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 20)
        self.assertEqual(actual[0], 'GGA')
        self.assertEqual(set(actual[1:]), set(('AA', 'aa', 'T', 'GGC')))
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 1)
        self.assertEqual(actual[0], 'GGA')
        self.assertEqual(set(actual[1:]), set(('AA', 'aa', 'T', 'GGC')))

    def test_same_as_ref_len(self):
        alt_seqs = ['AA', 'GGA', 'aa', 'GGA', 'T', 'GGC', 'aa']
        actual = viral_ngs.assemble.mummer.contig_chooser(alt_seqs, 1)
        self.assertEqual(actual[0], 'T')


class TestWgsimTool(TestCaseWithTmp):
    ''' Test wgsim tool wrapper helper methods '''

    def test_slice_fasta_whole_sequence(self):
        ''' Test slicing a fasta to a specific sequence ID without coordinate range '''
        wgsim = assemble.wgsim.WgsimTool()

        # Create test fasta with multiple sequences
        in_fasta = viral_ngs.core.file.mkstempfname('.fasta')
        seqs = [('seq1', 'ACGTACGTACGT'), ('seq2', 'GGGGCCCCAAAA'), ('seq3', 'TTTTTTTTTTTT')]
        makeFasta(seqs, in_fasta)

        # Slice to seq2
        out_fasta = viral_ngs.core.file.mkstempfname('.sliced.fasta')
        wgsim.slice_fasta(in_fasta, out_fasta, seq_id='seq2')

        # Verify only seq2 is in output
        records = list(Bio.SeqIO.parse(out_fasta, 'fasta'))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, 'seq2')
        self.assertEqual(str(records[0].seq), 'GGGGCCCCAAAA')

    def test_slice_fasta_with_coordinates(self):
        ''' Test slicing a fasta to a specific coordinate range '''
        wgsim = assemble.wgsim.WgsimTool()

        # Create test fasta
        in_fasta = viral_ngs.core.file.mkstempfname('.fasta')
        seqs = [('seq1', 'ACGTACGTACGTACGT')]
        makeFasta(seqs, in_fasta)

        # Slice to positions 5-10 (1-based inclusive)
        out_fasta = viral_ngs.core.file.mkstempfname('.sliced.fasta')
        wgsim.slice_fasta(in_fasta, out_fasta, seq_id='seq1', start=5, end=10)

        # Verify sliced sequence
        records = list(Bio.SeqIO.parse(out_fasta, 'fasta'))
        self.assertEqual(len(records), 1)
        self.assertEqual(str(records[0].seq), 'ACGTAC')

    def test_slice_fasta_no_params(self):
        ''' Test that slice_fasta with no optional params is a no-op '''
        wgsim = assemble.wgsim.WgsimTool()

        # Create test fasta
        in_fasta = viral_ngs.core.file.mkstempfname('.fasta')
        seqs = [('seq1', 'ACGTACGT'), ('seq2', 'GGGGCCCC')]
        makeFasta(seqs, in_fasta)

        # Slice with no params
        out_fasta = viral_ngs.core.file.mkstempfname('.sliced.fasta')
        wgsim.slice_fasta(in_fasta, out_fasta)

        # Verify all sequences are in output (no-op)
        records = list(Bio.SeqIO.parse(out_fasta, 'fasta'))
        self.assertEqual(len(records), 2)
        self.assertEqual(str(records[0].seq), 'ACGTACGT')
        self.assertEqual(str(records[1].seq), 'GGGGCCCC')

    def test_coverage_to_read_pairs(self):
        ''' Test coverage to read pairs calculation '''
        wgsim = assemble.wgsim.WgsimTool()

        # Test case 1: 20X coverage, 1000bp sequence, 150bp reads
        # Expected: (20 * 1000) / (2 * 150) = 66.67 -> 67 pairs
        pairs = wgsim.coverage_to_read_pairs(20, 1000, 150)
        self.assertEqual(pairs, 67)

        # Test case 2: 10X coverage, 5000bp sequence, 100bp reads
        # Expected: (10 * 5000) / (2 * 100) = 250 pairs
        pairs = wgsim.coverage_to_read_pairs(10, 5000, 100)
        self.assertEqual(pairs, 250)

        # Test case 3: Fractional coverage (0.5X)
        # Expected: (0.5 * 1000) / (2 * 150) = 1.67 -> 2 pairs
        pairs = wgsim.coverage_to_read_pairs(0.5, 1000, 150)
        self.assertEqual(pairs, 2)

        # Test case 4: High coverage
        # Expected: (100 * 3000) / (2 * 150) = 1000 pairs
        pairs = wgsim.coverage_to_read_pairs(100, 3000, 150)
        self.assertEqual(pairs, 1000)


class TestSimulateIlluminaReads(TestCaseWithTmp):
    ''' Test simulate_illumina_reads command for generating synthetic reads '''

    def setUp(self):
        super(TestSimulateIlluminaReads, self).setUp()
        self.in_fasta = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ref.lasv.fasta')
        # Calculate expected sequence lengths
        self.seq_lengths = {}
        for record in Bio.SeqIO.parse(self.in_fasta, 'fasta'):
            self.seq_lengths[record.id] = len(record.seq)

    def count_bam_reads(self, bam_file):
        ''' Count the number of reads in a BAM file '''
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        read_count = int(subprocess.check_output(
            [samtools.install_and_get_path(), 'view', '-c', bam_file]
        ).decode('UTF-8').strip())
        return read_count

    def expected_read_count(self, coverage, seq_length, read_length=150):
        ''' Calculate expected number of individual reads (not pairs) '''
        # Each read pair contributes 2 reads
        pairs = assemble.wgsim.WgsimTool().coverage_to_read_pairs(
            coverage, seq_length, read_length
        )
        return pairs * 2

    def test_simulate_uniform_coverage(self):
        ''' Test simulating reads with uniform coverage across all sequences '''
        out_bam = viral_ngs.core.file.mkstempfname('.bam')
        coverage = 20.0
        read_length = 150

        # Run simulate_illumina_reads with uniform coverage
        args = [self.in_fasta, out_bam, str(coverage), '--read_length', str(read_length), '--random_seed', '12345']
        args = viral_ngs.assembly.parser_simulate_illumina_reads(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # Verify BAM was created
        self.assertTrue(os.path.isfile(out_bam))
        self.assertGreater(os.path.getsize(out_bam), 0)

        # Count reads in output
        actual_reads = self.count_bam_reads(out_bam)

        # Calculate expected reads (sum across all sequences)
        expected_reads = sum(
            self.expected_read_count(coverage, seq_len, read_length)
            for seq_len in self.seq_lengths.values()
        )

        # Verify within 10% tolerance
        tolerance = 0.10
        self.assertGreater(actual_reads, expected_reads * (1 - tolerance))
        self.assertLess(actual_reads, expected_reads * (1 + tolerance))

    def test_simulate_per_sequence_coverage(self):
        ''' Test simulating reads with per-sequence coverage specification '''
        out_bam = viral_ngs.core.file.mkstempfname('.bam')
        read_length = 150

        # Get sequence IDs
        seq_ids = list(self.seq_lengths.keys())
        # Specify different coverage for each sequence
        coverage_spec = f"{seq_ids[0]}:20x {seq_ids[1]}:5x"

        # Run simulate_illumina_reads with per-sequence coverage
        args = [self.in_fasta, out_bam, coverage_spec, '--read_length', str(read_length), '--random_seed', '12345']
        args = viral_ngs.assembly.parser_simulate_illumina_reads(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # Verify BAM was created
        self.assertTrue(os.path.isfile(out_bam))
        self.assertGreater(os.path.getsize(out_bam), 0)

        # Count reads in output
        actual_reads = self.count_bam_reads(out_bam)

        # Calculate expected reads
        expected_reads = (
            self.expected_read_count(20, self.seq_lengths[seq_ids[0]], read_length) +
            self.expected_read_count(5, self.seq_lengths[seq_ids[1]], read_length)
        )

        # Verify within 10% tolerance
        tolerance = 0.10
        self.assertGreater(actual_reads, expected_reads * (1 - tolerance))
        self.assertLess(actual_reads, expected_reads * (1 + tolerance))

    def test_simulate_bed_coverage(self):
        ''' Test simulating reads with BED file coverage specification '''
        out_bam = viral_ngs.core.file.mkstempfname('.bam')
        bed_file = viral_ngs.core.file.mkstempfname('.bed')
        read_length = 150

        # Get sequence IDs
        seq_ids = list(self.seq_lengths.keys())

        # Create BED file with different coverage for different regions
        with open(bed_file, 'w') as f:
            # First sequence: 10x coverage for positions 1-2000
            f.write(f"{seq_ids[0]}\t0\t2000\t.\t10\n")
            # Second sequence: 25x coverage for positions 1-3000
            f.write(f"{seq_ids[1]}\t0\t3000\t.\t25\n")

        # Run simulate_illumina_reads with BED file
        args = [self.in_fasta, out_bam, bed_file, '--read_length', str(read_length), '--random_seed', '12345']
        args = viral_ngs.assembly.parser_simulate_illumina_reads(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # Verify BAM was created
        self.assertTrue(os.path.isfile(out_bam))
        self.assertGreater(os.path.getsize(out_bam), 0)

        # Count reads in output
        actual_reads = self.count_bam_reads(out_bam)

        # Calculate expected reads
        expected_reads = (
            self.expected_read_count(10, 2000, read_length) +
            self.expected_read_count(25, 3000, read_length)
        )

        # Verify within 10% tolerance
        tolerance = 0.10
        self.assertGreater(actual_reads, expected_reads * (1 - tolerance))
        self.assertLess(actual_reads, expected_reads * (1 + tolerance))


class TestFastaTrimTerminalAmbigs(TestCaseWithTmp):
    ''' Test fasta-trim-terminal-ambigs.pl script functionality '''

    def test_script_runs_successfully(self):
        ''' Test that fasta-trim-terminal-ambigs.pl can run without errors '''
        # Create a test FASTA with terminal ambiguous bases
        in_fasta = viral_ngs.core.file.mkstempfname('.fasta')
        out_fasta = viral_ngs.core.file.mkstempfname('.trimmed.fasta')

        # Create sequences with terminal Ns and ambiguous bases
        test_seqs = [
            ('seq1', 'NNNNNATCGATCGATCGNNNN'),
            ('seq2', 'WWWWATCGATCGATCGYYYYY'),
            ('seq3', 'RRRRATCGATCGATCGSSS'),
        ]
        makeFasta(test_seqs, in_fasta)

        # Run fasta-trim-terminal-ambigs.pl
        script_path = os.path.join(os.path.dirname(__file__), '..', '..', 'fasta-trim-terminal-ambigs.pl')
        cmd = [script_path, in_fasta]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())

        # Check that script runs successfully
        self.assertEqual(result.returncode, 0,
                        f"Script failed with stderr: {result.stderr}\nstdout: {result.stdout}")

        # Verify output contains trimmed sequences
        self.assertIn('ATCGATCGATCG', result.stdout)

    def test_script_with_3rules_option(self):
        ''' Test fasta-trim-terminal-ambigs.pl with --3rules option '''
        # Use the ebov-makona.fasta from test inputs
        in_fasta = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta')

        # Run with --3rules and standard parameters
        script_path = os.path.join(os.path.dirname(__file__), '..', '..', 'fasta-trim-terminal-ambigs.pl')
        cmd = [
            script_path,
            '--3rules',
            '--ten', '4',
            '--fifty', '15',
            '--maxfrac', '0.9',
            in_fasta
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())

        # Check that script runs successfully
        self.assertEqual(result.returncode, 0,
                        f"Script failed with stderr: {result.stderr}\nstdout: {result.stdout}")

        # Verify output was produced
        self.assertGreater(len(result.stdout), 0, "No output produced")
        self.assertIn('>', result.stdout, "No FASTA header in output")
