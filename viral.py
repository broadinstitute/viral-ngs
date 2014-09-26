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
import util.cmd, util.files, util.vcf

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
		for line in util.files.fastaMaker([(header, outseq)]):
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
		for header, seq in vcf_to_seqs(read_tabfile(args.inVcf),
			chrlens, samples, min_dp=args.min_dp, major_cutoff=args.major_cutoff,
			min_dp_ratio=args.min_dp_ratio):
			if args.trim_ends:
				seq = seq.strip('Nn')
			for line in util.files.fastaMaker([(header, seq)]):
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
				for line in util.files.fastaMaker([(record.id, str(record.seq).strip('Nn'))]):
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
				for line in util.files.fastaMaker([(record.id, ''.join(map(deambig_base, str(record.seq))))]):
					outf.write(line)
	log.info("done")
	return 0
__commands__.append(('deambig_fasta', main_deambig_fasta, parser_deambig_fasta))


def vcf_dpdiff(vcfs):
	for vcf in vcfs:
		samples = util.vcf.vcf_sample_names(vcf)
		assert len(samples)==1
		for row in read_tabfile(vcf):
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



def tabfile_values_rename(inFile, mapFile, outFile, col=0):
	# read map
	with open(mapFile, 'rt') as inf:
		name_map = dict(line.strip().split('\t') for line in inf)
	# convert file
	with open(outFile, 'wt') as outf:
		with open(inFile, 'rt') as inf:
			# copy header row verbatim
			outf.write(inf.readline())
			# all other rows: remap the specified column's values
			for line in inf:
				row = line.rstrip('\n').split('\t')
				row[col] = name_map[row[col]]
				outf.write('\t'.join(row)+'\n')
def parser_tabfile_rename():
	parser = argparse.ArgumentParser(
		description='''Take input tab file and copy to an output file while changing
the values in a specific column based on a mapping file.  The first line will pass
through untouched (it is assumed to be a header).''')
	parser.add_argument("inFile", help="Input flat file")
	parser.add_argument("mapFile",
		help="""Map file.  Two-column headerless file that maps input values to
		output values.  This script will error if there are values in inFile that do
		not exist in mapFile.""")
	parser.add_argument("outFile", help="Output flat file")
	parser.add_argument("--col_idx", dest="col_idx", type=int,
		help="""Which column number to replace (0-based index). [default: %(default)s]""",
		default=0)
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_tabfile_rename(args):
	tabfile_values_rename(args.inFile, args.mapFile, args.outFile, args.col_idx)
	return 0
__commands__.append(('tabfile_rename', main_tabfile_rename, parser_tabfile_rename))
			


def pos_to_number(row):
	row['pos'] = int(float(row['pos']))
	return row
def reposition_vphaser_deletions(row):
	if row['var'].startswith('D'):
		for k in ('ct_1','ct_2','ct_3','ct_4','extra1','extra2'):
			assert row.get(k,'D')[0] in ('D','i')
		row['pos'] = row['pos']-1
	return row

def vphaser_to_vcf(inFile, refFasta, multiAlignment, outVcf):
	''' Convert vPhaser2 parsed filtered output text file into VCF format '''
	
	# read in multiple alignments of consensus sequences
	aln = Bio.AlignIO.read(open(multiAlignment, 'rt'), 'fasta')

	# open reference genome and set ref as a BioPython SeqRecord
	ref = list(Bio.SeqIO.parse(open(refFasta, 'rt'), 'fasta'))
	assert len(ref)==1
	ref = ref[0]
	
	# prepare sample list
	samples = list(unique(row['patient'] for row in read_tabfile_dict(inFile)))
	samples_assembled = [(i, seq.id.split('.')[0], seq.id) for i,seq in enumerate(aln)]
	sample_idx_map = {}
	for s in samples:
		idx = [i for i,s_root,s_full in samples_assembled if s_root==s.split('.')[0]]
		assert len(idx)==1, "unable to uniquely find %s in %s" % (s, multiAlignment)
		sample_idx_map[s] = idx[0]
	
	# write output VCF file
	with open(outVcf, 'wt') as outf:
		outf.write('##fileformat=VCFv4.1\n')
		outf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		outf.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
		outf.write('##contig=<ID=%s,length=%d>\n' % (ref.id, len(ref)))
		outf.write('##reference=file://%s\n' % refFasta)
		header = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples
		outf.write('#'+'\t'.join(header)+'\n')
		
		# read in iSNVs and group rows based on unique position
		
		data = sorted(map(reposition_vphaser_deletions, map(pos_to_number, read_tabfile_dict(inFile))), key=lambda row: row['pos'])
		for pos, rows in itertools.groupby(data, lambda row: row['pos']):
			# get the set of alleles seen per patient
			rows = [(row['patient'], [row[h].split(':') for h in ('ct_1','ct_2','ct_3','ct_4','extra1','extra2') if row.get(h)]) for row in rows]
			# convert (allele, forward count, reverse count) tuples from strings to ints
			rows = [(s,[(a,int(f),int(r)) for a,f,r in counts]) for s,counts in rows]
			
			# filter based on same criteria used earlier (actually remove these calls)
			rows = [(s, list([(a,f,r) for a,f,r in counts if f>=5 and r>=5 and 10>=(float(f)/r)>=0.1]))
				for s,counts in rows]

			# remove patients where no internal variation exists anymore
			# remove (skip) positions where no patients exist anymore
			dropped = set(s for s,counts in rows if len(counts)<=1)
			rows = [(s,counts) for s,counts in rows if s not in dropped]
			if not rows:
				log.warn("dropping position %d due to loss of all samples" % pos)
				continue
			if dropped:
				log.warn("dropping samples %s at position %d due to filtered variation" % (dropped, pos))
			
			# combine fwd+rev counts and sort (allele,count) tuples in descending count order
			rows = [(s,list(sorted([(a,f+r) for a,f,r in counts], key=lambda (a,n):n, reverse=True))) for s,counts in rows]
			
			# define the length of this variation based on the largest deletion
			end = pos
			for s,counts in rows:
				for a,n in counts:
					if a.startswith('D'):
						end = max(end, pos+int(a[1:]))
			
			# find reference allele and consensus alleles
			refAllele = str(ref[pos-1:end].seq)
			consAlleles = dict((s, str(aln[sample_idx_map[s]][pos-1:end].seq)) for s in samples)
			for s,allele in consAlleles.items():
				if [a for a in allele if a not in set(('A','C','T','G'))]:
					log.warn("dropping unclean consensus for %s at %s-%s: %s" % (s, pos, end, allele))
					del consAlleles[s]
			
			# define genotypes and fractions
			iSNVs = {}
			rows = dict(rows)
			for s in samples:
				if s in rows:
					consAllele = consAlleles[s]
					# we have iSNV data on this sample
					tot_n = sum(n for a,n in rows[s])
					iSNVs[s] = {}
					for a,n in rows[s]:
						f = float(n)/tot_n
						if a.startswith('I'):
							# insertion allele is first ref base, plus inserted bases, plus subsequent ref bases
							a = consAllele[0] + a[1:] + consAllele[1:]
						elif a.startswith('D'):
							# deletion is the first ref base, plus remaining ref seq with the first few positions dropped off
							a = consAllele[0] + consAllele[1+int(a[1:]):]
						elif a in ('i','d'):
							# this is vphaser's way of saying the "reference" (majority/consensus) allele, in the face of other indel variants
							a = consAllele
						else:
							# this is a SNP
							assert a in set(('A','C','T','G'))
							if f>0.5 and a!=consAllele[0]:
								log.warn("vPhaser and assembly pipelines mismatch at %d/%s - consensus %s, vPhaser %s" % (pos, s, consAllele[0], a))
							a = a + consAllele[1:]
						assert a and a==a.upper()
						iSNVs[s][a] = f
					if unique(map(len, iSNVs[s].keys())) == [1]:
						assert consAllele in iSNVs[s].keys()
				elif s in consAlleles:
					# there is no iSNV data for this sample, so substitute the consensus allele
					iSNVs[s] = {consAlleles[s]:1.0}
			
			# get unique allele list and map to numeric
			alleles = [a for a,n in sorted(util.files.histogram(consAlleles.values()).items(), key=lambda(a,n):n, reverse=True) if a!=refAllele]
			alleles2 = list(itertools.chain(*[iSNVs[s].keys() for s in samples if s in iSNVs]))
			alleles = list(unique([refAllele] + alleles + alleles2))
			assert len(alleles)>1
			alleleMap = dict((a,i) for i,a in enumerate(alleles))
			genos = [str(alleleMap.get(consAlleles.get(s),'.')) for s in samples]
			freqs = [(s in iSNVs) and ','.join(map(str, [iSNVs[s].get(a,0.0) for a in alleles[1:]])) or '.' for s in samples]
			
			# prepare output row and write to file
			out = [ref.id, pos, '.', alleles[0], ','.join(alleles[1:]), '.', '.', '.', 'GT:AF']
			out = out + list(map(':'.join, zip(genos, freqs)))
			outf.write('\t'.join(map(str, out))+'\n')
			

def parser_vphaser_to_vcf():
	parser = argparse.ArgumentParser(
		description='''Convert vPhaser2 parsed filtered output text file into VCF format.
		We require the consensus assemblies for all these samples in a multi-alignment
		FASTA format as well, in order to resolve the ambiguity in vPhaser's output.
		All sample names and coordinates must be identical between inFile, inRef, and
		multiAlign.  We also require the reference genome FASTA (inRef) to determine
		reference alleles.  Requires a single-chromosome genome.''')
	parser.add_argument("inFile", help="Input vPhaser2 text file")
	parser.add_argument("inRef", help="Reference genome FASTA")
	parser.add_argument("multiAlign", help="Consensus genomes multi-alignment FASTA")
	parser.add_argument("outVcf", help="Output VCF file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_vphaser_to_vcf(args):
	vphaser_to_vcf(args.inFile, args.inRef, args.multiAlign, args.outVcf)
	return 0
__commands__.append(('vphaser_to_vcf', main_vphaser_to_vcf, parser_vphaser_to_vcf))



def compute_Fws(vcfrow):
	format = vcfrow[8].split(':')
	if 'AF' not in format:
		return None
	af_idx = format.index('AF')
	
	freqs = [dat.split(':') for dat in vcfrow[9:]]
	freqs = [float(dat[af_idx].split(',')[0]) for dat in freqs if len(dat)>af_idx and dat[af_idx]!='.' and dat[0]!='.' and int(dat[0])<=1]
	
	if len(freqs)<2:
		return None
	
	p_s = sum(freqs)/len(freqs)
	H_s = 2 * p_s * (1.0-p_s)
	
	if H_s==0.0:
		return None
	
	H_w = [2*p*(1.0-p) for p in freqs]
	H_w = sum(H_w)/len(H_w)
	return (H_s, 1.0 - H_w / H_s)

def add_Fws_vcf(inVcf, outVcf):
	with open(outVcf, 'wt') as outf:
		with util.files.open_or_gzopen(inVcf, 'rt') as inf:
			for line in inf:
				if line.startswith('##'):
					outf.write(line)
				elif line.startswith('#'):
					outf.write('##INFO=<ID=PI,Number=1,Type=Float,Description="Heterozygosity for this SNP in this sample set">\n')
					outf.write('##INFO=<ID=FWS,Number=1,Type=Float,Description="Fws statistic for iSNV to SNP comparisons (Manske 2012, Nature)">\n')
					outf.write(line)
				else:
					row = line.strip('\n').split('\t')
					Fws = compute_Fws(row)
					if Fws!=None:
						row[7] = row[7] + ";PI=%s;FWS=%s" % Fws
					outf.write('\t'.join(row)+'\n')

def parser_Fws():
	parser = argparse.ArgumentParser(
		description='''Compute the Fws statistic on iSNV data.
		See Manske, 2012 (Nature)''')
	parser.add_argument("inVcf", help="Input VCF file")
	parser.add_argument("outVcf", help="Output VCF file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_Fws(args):
	add_Fws_vcf(args.inVcf, args.outVcf)
	return 0
__commands__.append(('Fws', main_Fws, parser_Fws))


def iSNV_table(vcf_iter):
	for row in vcf_iter:
		info = dict(kv.split('=') for kv in row['INFO'].split(';'))
		samples = [k for k in row.keys() if k not in set(('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'))]
		for s in samples:
			f = row[s].split(':')[1]
			if f and f!='.':
				freqs = list(map(float, f.split(',')))
				f = sum(freqs)
				Hw = 1.0 - sum(p*p for p in [1.0-f]+freqs)
				out = {'chr':row['CHROM'], 'pos':row['POS'],
					'alleles':"%s,%s" %(row['REF'],row['ALT']), 'sample':s,
					'iSNV_freq':f, 'Hw':Hw}
				if 'EFF' in info:
					effs = [eff.rstrip(')').replace('(','|').split('|') for eff in info['EFF'].split(',')]
					effs = [[eff[i] for i in (0,3,4,5,6,9,11)] for eff in effs]
					effs = [eff for eff in effs if eff[5] not in ('sGP','ssGP') and int(eff[6])<2]
					assert len(effs)==1, "error at %s: %s" % (out['pos'], str(effs))
					eff = effs[0]
					if eff[2]:
						aa = eff[2].split('/')[0]
						assert aa.startswith('p.')
						aa = aa[2:]
						m = re.search(r"(\d+)", aa)
						out['eff_aa_pos'] = int(m.group(1))
					(out['eff_type'], out['eff_codon_dna'], out['eff_aa'], out['eff_prot_len'], out['eff_gene'], out['eff_protein'], rank) = eff
				if 'PI' in info:
					out['Hs_snp'] = info['PI']
				if 'FWS' in info:
					out['Fws_snp'] = info['FWS']
				yield out

def parser_iSNV_table():
	parser = argparse.ArgumentParser(
		description='''Convert VCF iSNV data to tabular text''')
	parser.add_argument("inVcf", help="Input VCF file")
	parser.add_argument("outFile", help="Output text file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_iSNV_table(args):
	header = ['pos','sample','patient','time','alleles','iSNV_freq','Hw',
		'eff_type','eff_codon_dna','eff_aa','eff_aa_pos','eff_prot_len','eff_gene','eff_protein']
	with open(args.outFile, 'wt') as outf:
		outf.write('\t'.join(header)+'\n')
		for row in iSNV_table(read_tabfile_dict(args.inVcf)):
			sample_parts = row['sample'].split('.')
			row['patient'] = sample_parts[0]
			if len(sample_parts)>1:
				row['time'] = sample_parts[1]
			outf.write('\t'.join(map(str, [row.get(h,'') for h in header]))+'\n')
	return 0
__commands__.append(('iSNV_table', main_iSNV_table, parser_iSNV_table))


def iSNP_per_patient(table, agg_fun=numpy.median):
	data = sorted(table, key=lambda row: (int(row['pos']), row['patient']))
	data = itertools.groupby(data, lambda row: (int(row['pos']), row['patient']))
	for x, rows in data:
		rows = list(rows)
		row = rows[0]
		if set(r['time'] for r in rows if r.get('time')):
			f = agg_fun(list(float(r['iSNV_freq']) for r in rows))
			row['iSNV_freq'] = f
			row['Hw'] = 2 * f * (1.0-f)
			row['sample'] = row['patient']
		else:
			assert len(rows)==1, "error, found multiple rows for %s:%s" % (row['pos'],row['patient'])
		yield row
def parser_iSNP_per_patient():
	parser = argparse.ArgumentParser(
		description='''Aggregate tabular iSNP data per patient x position (all time points averaged)''')
	parser.add_argument("inFile", help="Input text file")
	parser.add_argument("outFile", help="Output text file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_iSNP_per_patient(args):
	header = ['pos','patient','alleles','iSNV_freq','Hw',
		'eff_type','eff_codon_dna','eff_aa','eff_aa_pos','eff_prot_len','eff_gene','eff_protein']
	with open(args.outFile, 'wt') as outf:
		outf.write('\t'.join(header)+'\n')
		for row in iSNP_per_patient(read_tabfile_dict(args.inFile)):
			outf.write('\t'.join(map(str, [row.get(h,'') for h in header]))+'\n')
	return 0
__commands__.append(('iSNP_per_patient', main_iSNP_per_patient, parser_iSNP_per_patient))


def read_tabfile_dict(inFile):
	''' Read a tab text file (possibly gzipped) and return contents as an iterator of dicts. '''
	with util.files.open_or_gzopen(inFile, 'rt') as inf:
		header = None
		for line in inf:
			row = line.rstrip('\n').split('\t')
			if line.startswith('#'):
				row[0] = row[0][1:]
				header = row
			elif header==None:
				header = row
			else:
				assert len(header)==len(row)
				yield dict((k,v) for k,v in zip(header, row) if v)

def read_tabfile(inFile):
	''' Read a tab text file (possibly gzipped) and return contents as an iterator of arrays. '''
	with util.files.open_or_gzopen(inFile, 'rt') as inf:
		for line in inf:
			if not line.startswith('#'):
				yield line.rstrip('\n').split('\t')

def unique(items):
	''' Return unique items in the same order as seen in the input. '''
	seen = set()
	for i in items:
		if i not in seen:
			seen.add(i)
			yield i

if __name__ == '__main__':
	util.cmd.main_argparse(__commands__, __doc__)
