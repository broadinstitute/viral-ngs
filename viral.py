#!/usr/bin/env python
'''This script contains a number of utilities for a viral analysis pipeline
for Lassa virus sequence analysis, primarily used by Kristian Andersen.

Requires python >= 2.7 and BioPython.  On the Broad cluster, it is known
to work with the Python-2.7 and Python-3.4 dotkits.'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import os, tempfile, sys, shutil, argparse, logging, random, itertools, re, numpy
import Bio.AlignIO, Bio.SeqIO, Bio.Data.IUPACData
import util.cmd, util.file, util.vcf

log = logging.getLogger(__name__)
global_tool_paths = {}


def parser_modify_contig():
	parser = argparse.ArgumentParser(
		description='''Modifies an input contig. Depending on the options
		selected, can replace N calls with reference calls, replace ambiguous
		calls with reference calls, trim to the length of the reference, replace
		contig ends with reference calls, and trim leading and trailing Ns.
		Author: rsealfon.''')
	parser.add_argument("input", help="input alignment of reference and contig (should contain exactly 2 sequences)")
	parser.add_argument("output", help="Destination file for modified contigs")
	parser.add_argument("ref", help="reference sequence name (exact match required)")
	parser.add_argument("-n", "--name",
		help="fasta header output name (default:existing header)",
		default=None)
	parser.add_argument("-cn", "--call-reference-ns",
		help="""should the reference sequence be called if there is an
		N in the contig and a more specific base in the reference (default: %(default)s)""",
		default=False, action="store_true", dest="call_reference_ns")
	parser.add_argument("-t", "--trim-ends",
		help="should ends of contig.fasta be trimmed to length of reference (default: %(default)s)",
		default=False, action="store_true", dest="trim_ends")
	parser.add_argument("-r5", "--replace-5ends",
		help="should the 5'-end of contig.fasta be replaced by reference (default: %(default)s)",
		default=False, action="store_true", dest="replace_5ends")
	parser.add_argument("-r3", "--replace-3ends",
		help="should the 3'-end of contig.fasta be replaced by reference (default: %(default)s)",
		default=False, action="store_true", dest="replace_3ends")
	parser.add_argument("-l", "--replace-length",
		help="length of ends to be replaced (if replace-ends is yes) (default: %(default)s)",
		default=10, type=int)
	parser.add_argument("-f", "--format",
		help="Format for input alignment (default: %(default)s)",
		default="fasta")
	parser.add_argument("-r", "--replace-end-gaps",
		help="Replace gaps at the beginning and end of the sequence with reference sequence (default: %(default)s)",
		default=False, action="store_true", dest="replace_end_gaps")
	parser.add_argument("-rn", "--remove-end-ns",
		help="Remove leading and trailing N's in the contig (default: %(default)s)",
		default=False, action="store_true", dest="remove_end_ns")
	parser.add_argument("-ca", "--call-reference-ambiguous",
		help="""should the reference sequence be called if the contig seq is ambiguous and
		the reference sequence is more informative & consistant with the ambiguous base
		(ie Y->C) (default: %(default)s)""",
		default=False, action="store_true", dest="call_reference_ambiguous")
	util.cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',None)))
	return parser
def main_modify_contig(args):
    aln = Bio.AlignIO.read(args.input, args.format)
    ref_idx = find_ref_idx(aln,args.ref)
    assert ref_idx >= 0, "reference name '%s' not in alignment" % args.ref
    assert ref_idx <= 1, "alignment contains more than 2 sequences"

    consensus_idx = (ref_idx + 1) % 2

    ref = list(str(aln[ref_idx].seq))
    consensus = list(str(aln[consensus_idx].seq))

    assert len(ref) == len(consensus), "improper alignment"

    if (args.name == None):
        args.name = aln[consensus_idx].name
    
    if args.remove_end_ns:
        consensus = do_remove_end_ns(consensus)
    if args.call_reference_ns:
        consensus = do_call_reference_ns(ref, consensus)
    if args.call_reference_ambiguous:
        consensus = do_call_reference_ambiguous(ref, consensus)
    if args.trim_ends:
        consensus = do_trim_ends(ref, consensus)
    if args.replace_end_gaps:
        consensus = do_replace_end_gaps(ref, consensus)
    if args.replace_5ends:
        consensus = do_replace_5ends(ref, consensus, args.replace_length)
    if args.replace_3ends:
        consensus = do_replace_3ends(ref, consensus, args.replace_length)
        
    print_output(args.output, args.name, consensus)
    return 0
__commands__.append(('modify_contig', main_modify_contig, parser_modify_contig))


def find_ref_idx(aln,ref):
    if (aln[0].name == ref):
        return 0
    elif (aln[1].name == ref):
        return 1
    return -1

def do_call_reference_ns(ref, consensus):
    print("populating N's from reference...")
    for i in range(len(ref)):
        if (consensus[i].upper() == "N"):
            consensus[i] = ref[i]
    return consensus

def do_call_reference_ambiguous(ref, consensus):
    print("populating ambiguous bases from reference...")
    for i in range(len(ref)):
        if (consensus[i].upper() == "K"):
            if (ref[i].upper() == "G" or ref[i].upper() == "T"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "M"):
            if (ref[i].upper() == "A" or ref[i].upper() == "C"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "R"):
            if (ref[i].upper() == "A" or ref[i].upper() == "G"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "Y"):
            if (ref[i].upper() == "C" or ref[i].upper() == "T"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "S"):
            if (ref[i].upper() == "C" or ref[i].upper() == "G"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "W"):
            if (ref[i].upper() == "A" or ref[i].upper() == "T"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "B"):
            if (ref[i].upper() == "C" or ref[i].upper() == "G" or ref[i].upper() == "T"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "V"):
            if (ref[i].upper() == "A" or ref[i].upper() == "C" or ref[i].upper() == "G"):
                consensus[i] = ref[i]
        if (consensus[i].upper() == "D"):
            if (ref[i].upper() == "A" or ref[i].upper() == "G" or ref[i].upper() == "T"):
                consensus[i] = ref[i]
    return consensus

def do_trim_ends(ref, consensus):
    print("trimming ends...")
    for i in range(len(ref)):
        if (ref[i] != "-"):
            break
        else:
            consensus[i] = "-"
    for i in range(len(ref)):
        if (ref[len(ref) - i - 1] != "-"):
            break
        else:
            consensus[i] = "-"
    return consensus

def do_replace_end_gaps(ref, consensus):
    print("populating leading and trailing gaps from reference...")
    for i in range(len(ref)):
        if (consensus[i] != "-"):
            break
        consensus[i] = ref[i]
    for i in range(len(ref)):
        if (consensus[len(ref) - i - 1] != "-"):
            break
        consensus[len(ref) - i - 1] = ref[len(ref) - i - 1]
    return consensus

def do_replace_5ends(ref, consensus, replace_length):
    print("replacing 5' ends...")
    ct = 0
    for i in range(len(ref)):
        if (ref[i] != "-"):
            ct = ct + 1
        if (ct == replace_length):
            for j in range(0, i+1):
                consensus[j] = ref[j]
            break
    return consensus

def do_replace_3ends(ref, consensus, replace_length):
    print("replacing 3' ends...")
    ct = 0
    for i in range(len(ref)):
        if (ref[len(ref) - i - 1] != "-"):
            ct = ct + 1
        if (ct == replace_length):
            for j in range(0, i+1):
                consensus[len(ref)-j-1] = ref[len(ref)-j-1]
            break
    return consensus

def do_remove_end_ns(consensus):
    print("removing leading and trailing N's...")
    for i in range(len(consensus)):
        if (consensus[i].upper() == "N" or consensus[i] == "-"):
            consensus[i] = "-"
        else:
            break
    for i in range(len(consensus)):
        if (consensus[len(consensus) - i - 1].upper() == "N" or consensus[len(consensus) - i - 1] == "-"):
            consensus[len(consensus) - i - 1] = "-"
        else:
            break
    return consensus

def print_output(outfile, header, consensus):
	with open(outfile, "wt") as f:
		outseq = [x for x in consensus if not '-' in x]
		for line in util.file.fastaMaker([(header, outseq)]):
			f.write(line)



class MutableSequence:
	def __init__(self, name, start, stop, init_seq=None):
		assert stop>=start>=1
		if init_seq==None:
			self.seq = list('N' * (stop-start+1))
		else:
			self.seq = list(init_seq)
		assert stop-start+1 == len(self.seq)
		self.start = start
		self.stop = stop
		self.name = name
	def modify(self, p, new_base):
		assert self.start <= p <= self.stop
		i = p-self.start
		self.seq[i] = new_base
	def emit(self):
		return (self.name, ''.join(self.seq))

def alleles_to_ambiguity(allelelist):
	''' Convert a list of DNA bases to a single ambiguity base.
		All alleles must be one base long.  '''
	for a in allelelist:
		assert len(a)==1
	if len(allelelist)==1:
		return allelelist[0]
	else:
		# if you want to eliminate the dependency on BioPython, use this hardcoded line instead of the following one:
		#convert = {('C', 'G', 'T'): 'B', ('T',): 'T', ('A',): 'A', ('A', 'T'): 'W', ('A', 'C', 'T'): 'H', ('C',): 'C', ('A', 'G', 'T'): 'D', ('A', 'C', 'G'): 'V', ('G', 'T'): 'K', ('A', 'G'): 'R', ('C', 'G'): 'S', ('G',): 'G', ('A', 'C', 'G', 'T'): 'N', ('A', 'C'): 'M', ('C', 'T'): 'Y'}
		convert = dict([(tuple(sorted(v)),k) for k,v in Bio.Data.IUPACData.ambiguous_dna_values.items() if k!='X'])
		key = tuple(sorted(set(a.upper() for a in allelelist)))
		return convert[key]

def vcfrow_parse_and_call_snps(vcfrow, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
	''' Parse a single row of a VCF file, emit an iterator over each sample,
		call SNP genotypes using custom viral method based on read counts.
	'''
	# process this row
	c = vcfrow[0]
	p = int(vcfrow[1])
	alleles = [vcfrow[3]] + [a for a in vcfrow[4].split(',') if a not in '.']
	format = vcfrow[8].split(':')
	format = dict((format[i], i) for i in range(len(format)))
	assert 'GT' in format and format['GT']==0  # required by VCF spec
	assert len(vcfrow)==9+len(samples)
	info = dict(x.split('=') for x in vcfrow[7].split(';') if x != '.')
	info_dp = int(info.get('DP',0))
	
	# process each sample
	for i in range(len(samples)):
		sample = samples[i]
		rec = vcfrow[i+9].split(':')
		# require a minimum read coverage
		if len(alleles)==1:
			# simple invariant case
			dp = ('DP' in format and len(rec)>format['DP']) and int(rec[format['DP']]) or 0
			if dp < min_dp:
				continue
			geno = alleles
			if info_dp and float(dp)/info_dp < min_dp_ratio:
				log.warn("dropping invariant call at %s:%s %s (%s) due to low DP ratio (%s / %s = %s < %s)" % (
					c,p,sample,geno,dp,info_dp,float(dp)/info_dp,min_dp_ratio))
				continue
		else:
			# variant: manually call the highest read count allele if it exceeds a threshold
			assert ('AD' in format and len(rec)>format['AD'])
			allele_depths = list(map(int, rec[format['AD']].split(',')))
			assert len(allele_depths)==len(alleles)
			allele_depths = [(allele_depths[i], alleles[i]) for i in range(len(alleles)) if allele_depths[i]>0]
			allele_depths = list(reversed(sorted((n,a) for n,a in allele_depths if n>=min_dp)))
			if not allele_depths:
				continue
			dp = sum(n for n,a in allele_depths)

			if allele_depths[0][0] > (dp*major_cutoff):
				# call a single allele at this position if it is a clear winner
				geno = [allele_depths[0][1]]
			else:
				# call multiple alleles at this position if there is no clear winner
				geno = [a for n,a in allele_depths]
		if geno:
			yield (c, p, sample, geno)

def vcf_to_seqs(vcfIter, chrlens, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
	''' Take a VCF iterator and produce an iterator of chromosome x sample full sequences.'''
	seqs = {}
	cur_c = None
	for vcfrow in vcfIter:
		try:
			for c,p,s,alleles in vcfrow_parse_and_call_snps(vcfrow, samples, min_dp=min_dp, major_cutoff=major_cutoff, min_dp_ratio=min_dp_ratio):
				# changing chromosome?
				if c != cur_c:
					if cur_c!=None:
						# dump the previous chromosome before starting a new one
						for s in samples:
							yield seqs[s].emit()
			
					# prepare base sequences for this chromosome
					cur_c = c
					for s in samples:
						name = len(samples)>1 and ("%s-%s" % (c,sample)) or c
						seqs[s] = MutableSequence(name, 1, chrlens[c])
		
				# modify sequence for this chromosome/sample/position
				seqs[s].modify(p, alleles_to_ambiguity(alleles))
		except:
			log.exception("Exception occurred while parsing VCF file.  Row: '%s'" % vcfrow)
			raise
	
	# at the end, dump the last chromosome
	if cur_c!=None:
		for s in samples:
			yield seqs[s].emit()


def parser_vcf_to_fasta():
	parser = argparse.ArgumentParser(
		description='''Take input genotypes (VCF) and construct a consensus sequence
		(fasta) by using majority-read-count alleles in the VCF.
		Genotypes in the VCF will be ignored--we will use the allele
		with majority read support (or an ambiguity base if there is no clear majority).
		Uncalled positions will be emitted as N's.  Author: dpark.''')
	parser.add_argument("inVcf", help="Input VCF file")
	parser.add_argument("outFasta", help="Output FASTA file")
	parser.add_argument("--trim_ends",
		action="store_true", dest="trim_ends",
		default=False,
		help="""If specified, we will strip off continuous runs of N's from the beginning
		and end of the sequences before writing to output.  Interior N's will not be
		changed.""")
	parser.add_argument("--min_coverage", dest="min_dp", type=int,
		help="""Specify minimum read coverage (with full agreement) to make a call.
		[default: %(default)s]""",
		default=3)
	parser.add_argument("--major_cutoff", dest="major_cutoff", type=float,
		help="""If the major allele is present at a frequency higher than this cutoff,
		we will call an unambiguous base at that position.  If it is equal to or below
		this cutoff, we will call an ambiguous base representing all possible alleles at
		that position. [default: %(default)s]""",
		default=0.5)
	parser.add_argument("--min_dp_ratio", dest="min_dp_ratio", type=float,
		help="""The input VCF file often reports two read depth values (DP)--one for
		the position as a whole, and one for the sample in question.  We can optionally
		reject calls in which the sample read count is below a specified fraction of the
		total read count.  This filter will not apply to any sites unless both DP values
		are reported.  [default: %(default)s]""",
		default=None)
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_vcf_to_fasta(args):
	assert args.min_dp >= 0
	assert 0.0 <= args.major_cutoff < 1.0
	
	chrlens = dict(util.vcf.vcf_chrlens(args.inVcf))
	samples = util.vcf.vcf_sample_names(args.inVcf)
	with open(args.outFasta, 'wt') as outf:
		for header, seq in vcf_to_seqs(util.file.read_tabfile(args.inVcf),
			chrlens, samples, min_dp=args.min_dp, major_cutoff=args.major_cutoff,
			min_dp_ratio=args.min_dp_ratio):
			if args.trim_ends:
				seq = seq.strip('Nn')
			for line in util.file.fastaMaker([(header, seq)]):
				outf.write(line)
	
	# done
	log.info("done")
	return 0
__commands__.append(('vcf_to_fasta', main_vcf_to_fasta, parser_vcf_to_fasta))



def parser_trim_fasta():
	parser = argparse.ArgumentParser(
		description='''Take input sequences (fasta) and trim any continuous sections of
		N's from the ends of them.  Write trimmed sequences to an output fasta file.''')
	parser.add_argument("inFasta", help="Input fasta file")
	parser.add_argument("outFasta", help="Output (trimmed) fasta file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_trim_fasta(args):
	with open(args.outFasta, 'wt') as outf:
		with open(args.inFasta, 'rt') as inf:
			for record in Bio.SeqIO.parse(inf, 'fasta'):
				for line in util.file.fastaMaker([(record.id, str(record.seq).strip('Nn'))]):
					outf.write(line)
	log.info("done")
	return 0
__commands__.append(('trim_fasta', main_trim_fasta, parser_trim_fasta))



def deambig_base(base):
	''' take a single base (possibly a IUPAC ambiguity code) and return a random
		non-ambiguous base from among the possibilities '''
	return random.choice(Bio.Data.IUPACData.ambiguous_dna_values[base.upper()])
	
def parser_deambig_fasta():
	parser = argparse.ArgumentParser(
		description='''Take input sequences (fasta) and replace any ambiguity bases with a
		random unambiguous base from among the possibilities described by the ambiguity
		code.  Write output to fasta file.''')
	parser.add_argument("inFasta", help="Input fasta file")
	parser.add_argument("outFasta", help="Output fasta file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_deambig_fasta(args):
	with open(args.outFasta, 'wt') as outf:
		with open(args.inFasta, 'rt') as inf:
			for record in Bio.SeqIO.parse(inf, 'fasta'):
				for line in util.file.fastaMaker([(record.id, ''.join(map(deambig_base, str(record.seq))))]):
					outf.write(line)
	log.info("done")
	return 0
__commands__.append(('deambig_fasta', main_deambig_fasta, parser_deambig_fasta))


def vcf_dpdiff(vcfs):
	for vcf in vcfs:
		samples = util.vcf.vcf_sample_names(vcf)
		assert len(samples)==1
		for row in util.file.read_tabfile(vcf):
			dp1 = int(dict(x.split('=') for x in row[7].split(';') if x != '.').get('DP',0))
			dp2 = 0
			if 'DP' in row[8].split(':'):
				dpidx = row[8].split(':').index('DP')
				if len(row[9].split(':'))>dpidx:
					dp2 = int(row[9].split(':')[dpidx])
			ratio = ''
			if dp1:
				ratio = float(dp2)/dp1
			yield (row[0],row[1],samples[0],dp1,dp2,dp1-dp2,ratio)
	

def parser_dpdiff():
	parser = argparse.ArgumentParser(
		description='''Take input VCF files (all with only one sample each) and report
		on the discrepancies between the two DP fields (one in INFO and one in the
		sample's genotype column).''')
	parser.add_argument("inVcfs", help="Input VCF file", nargs='+')
	parser.add_argument("outFile", help="Output flat file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_dpdiff(args):
	header = ['chr','pos','sample','dp_info','dp_sample','diff','ratio']
	with open(args.outFile, 'wt') as outf:
		outf.write('#'+'\t'.join(header)+'\n')
		for row in vcf_dpdiff(args.inVcfs):
			outf.write('\t'.join(map(str, row))+'\n')
	return 0
__commands__.append(('dpdiff', main_dpdiff, parser_dpdiff))


if __name__ == '__main__':
	util.cmd.main_argparse(__commands__, __doc__)
