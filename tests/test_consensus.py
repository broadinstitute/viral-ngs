# Unit tests for consensus.py

__author__ = "dpark@broadinstitute.org"

import consensus, util.cmd, util.file
import Bio.SeqIO, Bio.Data.IUPACData
import unittest
import os, shutil, tempfile, argparse, itertools


def set_tmpDir(name):
	proposed_prefix = ['tmp']
	if name:
		proposed_prefix.append(name)
	for e in ('LSB_JOBID','LSB_JOBINDEX'):
		if e in os.environ:
			proposed_prefix.append(os.environ[e])
	tempfile.tempdir = tempfile.mkdtemp(prefix='-'.join(proposed_prefix)+'-',
		dir=util.cmd.find_tmpDir())
def destroy_tmpDir():
	shutil.rmtree(tempfile.tempdir)

def makeFasta(seqs, outFasta):
	with open(outFasta, 'wt') as outf:
		for line in util.file.fastaMaker(seqs):
			outf.write(line)

class TestCommandHelp(unittest.TestCase):
	def test_help_parser_for_each_command(self):
		for cmd_name, main_fun, parser_fun in consensus.__commands__:
			parser = parser_fun()
			helpstring = parser.format_help()


class TestAmbiguityBases(unittest.TestCase):
	def test_non_failure(self):
		''' Make sure that alleles_to_ambiguity runs without errors for every possible
			combination of inputs.  Check that the output is one-character long and uppercase.
		'''
		bases = ('A','C','T','G')
		for i in range(1,5):
			for alleles in itertools.permutations(bases, i):
				out = consensus.alleles_to_ambiguity(alleles)
				self.assertEqual(1, len(out))
				self.assertEqual(out, out.upper())



class TestManualSnpCaller(unittest.TestCase):
	''' Test the vcfrow_parse_and_call_snps method.. lots of edge cases. '''
	def test_missing_dp(self):
		''' VCF files might contain rows with no calls or any kind of data and that's okay. '''
		row = ['chr10', '105', '.', 'G', '.', '.', '.', '.', 'GT', './.']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
		self.assertEqual(out, [])
	def test_dp_inaccurate(self):
		''' The DP might not equal the sum of the ADs and that's okay apparently. '''
		row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:5:2,2']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
		self.assertEqual(set(out[0][3]), set(['G','A']))
		row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:2:3,3']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(set(out[0][3]), set(['G','A']))
		row = ['chr10', '105', '.', 'G', 'A', '.', '.', '.', 'GT:DP:AD', '0/1/1:10:2,0']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(out, [])
	def test_invariant_sites(self):
		''' Invariant site handling is slightly different in code, so test it specially. '''
		row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', '0/0:3']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(out, [('LASV.l',1,'s1',['T'])])
		row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=0))
		self.assertEqual(out, [('LASV.l',1,'s1',['T'])])
		row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', '0/0']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
		self.assertEqual(out, [])
		row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT', './.']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
		self.assertEqual(out, [])
		row = ['LASV.l', '1', '.', 'T', '.', '.', '.', '.', 'GT:DP', './.:10']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=1))
		self.assertEqual(out, [('LASV.l',1,'s1',['T'])])
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
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(set(out[0][3]), set(['G','A','C']))
		row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,3,0,3']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(set(out[0][3]), set(['A','T']))
		row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(out, [])
		row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,0,2']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=2))
		self.assertEqual(set(out[0][3]), set(['A','T']))
		row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:2,0,3,0']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(out[0][3], ['C'])
		row = ['thecontig', '105000', '.', 'G', 'A,C,T', '.', '.', '.', 'GT:AD', '0/1:0,2,3,4']
		out = list(consensus.vcfrow_parse_and_call_snps(row, ['s1'], min_dp=3))
		self.assertEqual(out[0][3], ['T'])
		


class TestDeambigAndTrimFasta(unittest.TestCase):
	''' Test the deambig_fasta and trim_fasta commands. '''
	def setUp(self):
		set_tmpDir('TestDeambigAndTrimFasta')
	def tearDown(self):
		destroy_tmpDir()
	def run_method(self, inseqs, parser_fun, main_fun):
		fasta_in = util.file.mkstempfname()
		fasta_out = util.file.mkstempfname()
		makeFasta([(str(i), inseqs[i]) for i in range(len(inseqs))], fasta_in)
		args = parser_fun().parse_args([fasta_in, fasta_out])
		main_fun(args)
		return (fasta_in, fasta_out)
	def test_trim_fasta(self):
		''' Simple test of the trim_fasta command '''
		inseqs = ['NNnnNNnNaslkdfjasdkfNNNN','NNNnnN','NNN123','ATCG']
		expected = ['aslkdfjasdkf','','123','ATCG']
		expected = dict((str(i), expected[i]) for i in range(len(expected)))
		fasta_in, fasta_out = self.run_method(inseqs, consensus.parser_trim_fasta, consensus.main_trim_fasta)
		with open(fasta_out, 'rt') as fa:
			for record in Bio.SeqIO.parse(fa, 'fasta'):
				self.assertIn(record.id, expected)
				self.assertEqual(str(record.seq), expected[record.id])
	def test_deambig_fasta(self):
		''' Simple test of the deambig_fasta command '''
		table = [(k,v) for k,v in Bio.Data.IUPACData.ambiguous_dna_values.items() if k!='X']
		keys = [k for k,v in table]
		vals = [set(v) for k,v in table]
		keys = keys + [k.lower() for k in keys]
		vals = vals + vals
		inseq = ''.join(keys)
		fasta_in, fasta_out = self.run_method([inseq], consensus.parser_deambig_fasta, consensus.main_deambig_fasta)
		with open(fasta_out, 'rt') as fa:
			for rec in Bio.SeqIO.parse(fa, 'fasta'):
				self.assertEqual(rec.id, '0')
				outseq = str(rec.seq)
				for i in range(len(outseq)):
					self.assertIn(outseq[i], vals[i])
