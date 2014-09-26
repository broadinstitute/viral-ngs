#!/usr/bin/env python
'''This script contains a number of utilities for intrahost variant calling
and annotation for viral genomes.

Requires python >= 2.7 and BioPython.  On the Broad cluster, it is known
to work with the Python-2.7 and Python-3.4 dotkits.'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org, swohl@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import os, tempfile, sys, shutil, argparse, logging, random, itertools, re, numpy
import Bio.AlignIO, Bio.SeqIO, Bio.Data.IUPACData
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)
global_tool_paths = {}


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
	samples = list(util.misc.unique(row['patient'] for row in util.file.read_tabfile_dict(inFile)))
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
		
		data = sorted(map(reposition_vphaser_deletions, map(pos_to_number, util.file.read_tabfile_dict(inFile))), key=lambda row: row['pos'])
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
					if util.misc.unique(map(len, iSNVs[s].keys())) == [1]:
						assert consAllele in iSNVs[s].keys()
				elif s in consAlleles:
					# there is no iSNV data for this sample, so substitute the consensus allele
					iSNVs[s] = {consAlleles[s]:1.0}
			
			# get unique allele list and map to numeric
			alleles = [a for a,n in sorted(util.misc.histogram(consAlleles.values()).items(), key=lambda(a,n):n, reverse=True) if a!=refAllele]
			alleles2 = list(itertools.chain(*[iSNVs[s].keys() for s in samples if s in iSNVs]))
			alleles = list(util.misc.unique([refAllele] + alleles + alleles2))
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
		with util.file.open_or_gzopen(inVcf, 'rt') as inf:
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
		for row in iSNV_table(util.file.read_tabfile_dict(args.inVcf)):
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
		for row in iSNP_per_patient(util.file.read_tabfile_dict(args.inFile)):
			outf.write('\t'.join(map(str, [row.get(h,'') for h in header]))+'\n')
	return 0
__commands__.append(('iSNP_per_patient', main_iSNP_per_patient, parser_iSNP_per_patient))


if __name__ == '__main__':
	util.cmd.main_argparse(__commands__, __doc__)
