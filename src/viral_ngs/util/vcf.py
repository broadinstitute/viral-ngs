'''This gives a number of useful quick methods for dealing with VCF files.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import os, shutil, logging, itertools, sqlite3
import pysam
import util.file, util.misc

log = logging.getLogger(__name__)

def make_intervals(i, n, fasta, chr_prefix='', verbose=False):
	''' Divide a sorted genome into n equally sized parts and return the i'th
		part.  We will return a list of intervals: chr, start, stop.  It may
		contain multiple chromosomes in order to fill out a total length equal
		to the other parts.  Each part will be adjacent and non-overlapping with
		the next part. i must be a number from 1 to n.
	'''
	assert 1 <= i <= n
	
	# read genome dict file
	tot = 0
	chrlens = []
	for c,c_len in get_chrlens(fasta):
		if c.startswith(chr_prefix):
			chrlens.append((c,c_len,tot))
			tot += c_len
	
	# define our chunk by gpos:
	part_size = tot//n
	g_start = 1 + part_size * (i-1)
	g_stop = part_size * i
	if i==n:
		g_stop = tot
	
	# find the genomic intervals that correspond to our gpos window
	out = []
	for c, c_len, c_g_start in chrlens:
		c_g_stop = c_g_start + c_len
		c_g_start += 1
		if c_g_stop >= g_start and c_g_start <= g_stop:
			start = max(g_start, c_g_start) - c_g_start + 1
			stop  = min(g_stop,  c_g_stop)  - c_g_start + 1
			out.append((c, start, stop))
	
	if verbose:
		log.info("Dividing the %d bp genome into %d chunks of %d bp each.  The %dth chunk contains the following %d intervals: %s" % (
			tot, n, part_size, i, len(out), ', '.join(["%s:%d-%d"%x for x in out])))
	return out


def sliding_windows(fasta, width, offset, chr_prefix=''):
	''' Divide a genome into sliding windows and return as an iterator of
		intervals (chr, start, stop).  The windows are of fixed width
		(except maybe the last window on each chromosome) and may overlap
		(offset<width) or be discontinuous (offset>width).
	'''
	assert width>0 and offset>0
	for c,c_len in get_chrlens(fasta):
		if c.startswith(chr_prefix):
			start = 1
			while start <= c_len:
				stop = min(c_len, start + width - 1)
				yield (c, start, stop)
				start += offset


def vcf_bgzip_index(inVcf, outVcf, tabixPath=None, vcftoolsPath=None):
	assert (inVcf.endswith('.vcf') or inVcf.endswith('.vcf.gz')) and outVcf.endswith('.vcf.gz')
	''' Compress (bgzip) and index (tabix and/or vcftools) a VCF file.'''
	
	if inVcf.endswith('.vcf'):
		log.info("compressing output to %s" % outVcf)
		cmdline = "%s/bgzip -c %s > %s" % (tabixPath, inVcf, outVcf)
		assert not os.system(cmdline)
	elif inVcf==outVcf:
		log.info("leaving %s, already compressed" % inVcf)
	else:
		log.info("copying compressed vcf to %s" % outVcf)
		shutil.copy(inVcf, outVcf)
	
	log.info("indexing with tabix")
	cmdline = "%s/tabix %s -f -p vcf" % (tabixPath, outVcf)
	assert not os.system(cmdline)
	
	if vcftoolsPath:
		log.info("indexing with vcftools")
		tmpFile = util.file.mkstempfname(prefix='vcftools-log-', suffix='.vcf')
		cmdline = "%s/vcftools --gzvcf %s --out %s --force-index-write" % (vcftoolsPath, outVcf, tmpFile)
		assert not os.system(cmdline)
		os.unlink(tmpFile)
		os.unlink(tmpFile+'.log')
	return outVcf

class GenomePosition:
	''' Provide a mapping from chr:pos to genomic position.
		Read chromosome lengths and order from either a Picard/GATK-index for
		a FASTA file (a .dict file) or from a VCF header.
	'''
	def __init__(self, seqDb):
		self.gpos_map = {}
		self.clen_map = {}
		self.chrs = []
		totlen = 0
		for c,clen in get_chrlens(seqDb):
			self.chrs.append((c,clen))
			self.gpos_map[c] = totlen
			self.clen_map[c] = clen
			totlen += clen
		self.total = totlen
	def get_gpos(self, c, p):
		assert isinstance(p, int)
		assert c in self.gpos_map
		assert 1 <= p <= self.clen_map[c]
		return p + self.gpos_map[c]
	def get_chr_pos(self, gpos):
		assert isinstance(gpos, int)
		assert 1 <= gpos <= self.total
		totlen = 0
		for c,clen in self.chrs:
			if gpos <= totlen+clen:
				break
			totlen += clen
		return (c,gpos-totlen)

def get_chrlens(inFile):
	''' Read chromosome lengths and order from either a Picard/GATK-index for
		a FASTA file (a .dict file) or from "contig" rows in the VCF header.
	'''
	chrlens = []
	if hasattr(inFile, 'chrlens'):
		chrlens = inFile.chrlens()
	else:
		if inFile.endswith('.fasta'):
			inFile = inFile[:-len('.fasta')] + '.dict'
		elif inFile.endswith('.fa'):
			inFile = inFile[:-len('.fa')] + '.dict'
		if inFile.endswith('.dict'):
			with open(inFile, 'rt') as inf:
				for line in inf:
					row = line.rstrip('\n').split('\t')
					if row[0]=='@SQ':
						assert row[1].startswith('SN:') and row[2].startswith('LN:')
						c = row[1][3:]
						c_len = int(row[2][3:])
						chrlens.append((c,c_len))
		elif inFile.endswith('.vcf') or inFile.endswith('.vcf.gz'):
			with util.file.open_or_gzopen(inFile, 'rt') as inf:
				for line in inf:
					line = line.rstrip('\n')
					if line.startswith('##contig=<ID=') and line.endswith('>'):
						line = line[13:-1]
						c = line.split(',')[0]
						clen = int(line.split('=')[1])
						chrlens.append((c,clen))
					elif line.startswith('#CHROM'):
						break
		else:
			raise AssertionError("unrecognized file type %s" % inFile)
	assert chrlens, "no sequence data found in %s % inFile"
	return chrlens

def get_chroms(inVcf, tabixPath=None):
	''' Get a list of unique chromosomes for this genome, in sort order.
		(Use tabix to do it quickly)
	'''
	tmpFile = util.file.mkstempfname(prefix='chrnames-', suffix='.txt')
	cmdline = "%s/tabix -l %s > %s" % (tabixPath, inVcf, tmpFile)
	assert not os.system(cmdline)
	chroms = []
	with open(tmpFile, 'rt') as inf:
		for line in inf:
			chroms.append(line.rstrip('\r\n'))
	os.unlink(tmpFile)
	return chroms

def vcf_sample_names(inVcf):
	''' Return the list of sample names in a given VCF file (quickly). '''
	samples = None
	with util.file.open_or_gzopen(inVcf, 'rt') as inf:
		for line in inf:
			if line.startswith('#CHROM'):
				row = line.rstrip('\r\n').split('\t')
				samples = row[9:]
				break  # we must properly close the file.. not sure if "return" does that
	assert samples, "error: no header line!"
	return samples

def vcf_subset(inVcf, c, start_stop=None, outVcf=None, keepHeader=False, tabixPath=None):
	''' Pull just a piece of a VCF file into a new VCF file (create a temp
		file if outVcf is not specified).
	'''
	if outVcf==None:
		outVcf = util.file.mkstempfname(prefix='vcf_subset-%s-'%c, suffix='.vcf')
	assert inVcf.endswith('.vcf.gz') and outVcf.endswith('.vcf')
	cmdline = "%s/tabix" % tabixPath
	if keepHeader:
		cmdline += ' -h'
	cmdline += ' %s %s' % (inVcf, c)
	if start_stop:
		cmdline += ':%d-%d' % start_stop
	cmdline += ' > %s' % outVcf
	assert not os.system(cmdline)
	return outVcf

def vcf_subset_rows(inVcf, c=None, start_stop=None, tabixPath=None):
	''' Pull just a piece of a VCF file and return as an iterator of 4-tuples '''
	if c:
		sub_vcf = vcf_subset(inVcf,c,start_stop=start_stop,keepHeader=True,tabixPath=tabixPath)
		for x in vcf_rows(sub_vcf):
			yield x
		os.unlink(sub_vcf)
	else:
		for x in vcf_rows(inVcf):
			yield x

def vcf_rows(inVcf):
	''' Read a VCF file and return contents as an iterator.
		Return each row as a 4-tuple:
			1: [first 9 columns]
			2: [all genotype columns (10th col and onward)]
			3: {dict of sample name : genotype column}
			4: [list of sample names]
	'''
	with util.file.open_or_gzopen(inVcf, 'rt') as inf:
		samples = []
		for line in inf:
			if not line.startswith('#'):
				row = line.rstrip('\n').split('\t')
				yield (row[:9],
					row[9:],
					dict([(samples[i],row[9+i]) for i in range(len(samples))]),
					samples)
			elif line.startswith('#CHROM'):
				samples = line.rstrip('\n').split('\t')[9:]

def vcf_haploid_iterator(inVcf, sample_list=None,
	drop_indels=False, drop_monomorphic=False, drop_multiallelic=False,
	interval=None, tabixPath=None):
	'''	Read a VCF file and return contents as an iterator with parsed
		contents (assuming haploid genotypes).  Each row is returned as
		a 4-tuple:
			1: chr
			2: pos (int)
			3: list of allele strings (in order)
			4: list of genotypes as 2-tuples: (sample, int allele)
	'''
	if tabixPath or interval==None:
		c, start_stop = (None,None)
		if interval:
			iparts = interval.split(':')
			c = iparts[0]
			if len(iparts)>1:
				start_stop = map(int, iparts[1].split('-'))
		for info,genos,genoMap,rowsamp in vcf_subset_rows(inVcf, c, start_stop=start_stop, tabixPath=tabixPath):
			c,p = info[:2]
			p = int(p)
			alleles = [info[3]] + info[4].split(',')
			if drop_monomorphic and len(alleles)==1:
				continue
			if drop_indels and not all([len(a)==1 for a in alleles]):
				continue
			if sample_list!=None:
				rowsamp = sample_list
			genos = [(s,int(genoMap[s][0])) for s in rowsamp if genoMap[s][0]!='.']
			n_alleles = len(set([a for s,a in genos]))
			if drop_monomorphic and n_alleles<2 or drop_multiallelic and n_alleles>2:
				continue
			yield (c, p, alleles, genos)
	else:
		with VcfReader(inVcf) as vcf:
			for c,p,alleles,genos in vcf.get_range(region=interval, as_strings=False):
				if drop_monomorphic and len(alleles)==1:
					continue
				if drop_indels and not all([len(a)==1 for a in alleles]):
					continue
				if sample_list!=None:
					genos = dict(genos)
					genos = [(s,genos[s]) for s in sample_list if s in genos]
				n_alleles = len(set([a for s,a in genos]))
				if drop_monomorphic and n_alleles<2 or drop_multiallelic and n_alleles>2:
					continue
				yield (c, p, alleles, genos)

def calc_maf(genos, ancestral=None, ploidy=1):
	# get list of alleles
	if ploidy==1:
		alleles = genos
	else:
		alleles = []
		for g in genos:
			g = g.split('/')
			assert len(g)==ploidy
			alleles += g
	
	# count up
	out = {'n_tot':len(alleles)}
	acounts = util.misc.histogram(alleles)
	alist = sorted([(n,a) for a,n in acounts.items()])
	
	# daf
	if ancestral != None:
		out['a_ancestral'] = ancestral
		derived = list(sorted([a for a in acounts.keys() if a!=ancestral]))
		out['a_derived'] = ','.join(derived)
		out['dac'] = sum(acounts[a] for a in derived)
		out['daf'] = out['n_tot'] and float(out['dac'])/out['n_tot'] or None
	
	# maf
	if out['n_tot']:
		out['a_major'] = alist[-1][1]
		out['a_minor'] = ','.join([a for n,a in alist[:-1]])
		out['mac'] = out['n_tot'] - alist[-1][0]
		out['maf'] = float(out['mac'])/out['n_tot']
	else:
		out['a_major'] = None
		out['a_minor'] = None
		out['mac'] = None
		out['maf'] = None
		
	return out


class SnpDb:
	''' Present genotype data as a sqlite database.  Load from a VCF file. '''
	def __init__(self, inVcf,
		enforce_unique_pos=True, sample_list=None,
		drop_indels=False, drop_monomorphic=False, drop_multiallelic=False,
		interval=None):
		assert inVcf.endswith('.vcf.gz')
		clens = get_chrlens(inVcf)
		self.clens = dict(clens)
		self.contigs = [c for c,l in clens]
		self.sample_names = vcf_sample_names(inVcf)
		root_fname = inVcf[:-7].split('/')[-1]
		self.dbFile = util.file.mkstempfname(prefix='%s-'%root_fname,suffix='.db')
		self.conn = sqlite3.connect(self.dbFile, isolation_level='DEFERRED')
		self.cur = self.conn.cursor()
		self.cur.execute("""create table snp (
			chr string not null,
			pos integer not null,
			alleles string not null)""")
		self.cur.execute("""create table geno (
			chr string not null,
			pos integer not null,
			sample string not null,
			allele integer not null)""")
		self.cur.execute("create %s index snp_idx on snp(chr,pos)" % (
			enforce_unique_pos and "unique" or ""))
		self.cur.execute("create %s index geno_idx on geno(chr,pos,sample)" % (
			enforce_unique_pos and "unique" or ""))
		self.conn.commit()
		imap = hasattr(itertools, 'imap') and itertools.imap or map  #py2 & py3 compatibility
		self.cur.executemany("insert into snp (chr,pos,alleles) values (?,?,?)",
			imap(lambda snp: (snp[0],snp[1],','.join(snp[2])),
				vcf_haploid_iterator(inVcf,
					sample_list=sample_list,
					drop_indels=drop_indels, drop_monomorphic=drop_monomorphic,
					drop_multiallelic=drop_multiallelic,
					interval=interval)))
		self.conn.commit()
		self.cur.executemany("insert into geno (chr,pos,sample,allele) values (?,?,?,?)",
			self._geno_iterator(vcf_haploid_iterator(inVcf, sample_list=sample_list,
				drop_indels=drop_indels, drop_monomorphic=drop_monomorphic,
				drop_multiallelic=drop_multiallelic,
				interval=interval)))
		self.conn.commit()
	def _geno_iterator(self, vcf_iterator):
		for c,p,alleles,genos in vcf_iterator:
			for s,a in genos:
				yield (c,p,s,a)
	def close(self):
		self.cur.close()
		self.conn.close()
		os.unlink(self.dbFile)
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def chroms(self):
		return self.contigs
	def samples(self):
		return self.sample_names
	def chrlens(self):
		return self.clens
	def get_snp_genos(self, c, p, as_strings=True):
		self.cur.execute("select sample,allele from geno where chr=? and pos=?", [c,p])
		out = [(s,a) for s,a in self.cur]
		if as_strings and out:
			self.cur.execute("select alleles from snp where chr=? and pos=?", [c,p])
			alleles = [x[0] for x in self.cur]
			if not len(alleles)==1:
				out = []
			else:
				alleles = alleles[0].split(',')
				out = [(s,alleles[a]) for s,a in out]
		return dict(out)
	def get_positions(self, c=None, start=None, stop=None):
		if start==None:
			start = 0
		if stop==None:
			stop = 1E15
		if c==None:
			self.cur.execute("select chr,pos from snp where chr=? and pos>=? and pos<=? order by chr,pos", [c,start,stop])
		else:
			self.cur.execute("select chr,pos from snp where chr=? and pos>=? and pos<=? order by chr,pos", [c,start,stop])
		return [(c,p,p) for c,p in self.cur]
	def get_range(self, c, start, stop, as_strings=True):
		'''	Query for a region and return an iterator of 4-tuples:
				1: chr
				2: pos (int)
				3: list of allele strings (in order)
				4: list of genotypes as 2-tuples: (sample, allele)
			If as_strings, the alleles will be actual alleles.  Otherwise,
			alleles will be integers.
		'''
		cur2 = self.conn.cursor()
		cur2.execute("select pos,alleles from snp where chr=? and pos>=? and pos<=? order by pos", [c,start,stop])
		for p,alleles in cur2:
			self.cur.execute("select sample,allele from geno where chr=? and pos=?", [c,p])
			genos = [(s,a) for s,a in self.cur]
			alleles = alleles.split(',')
			if as_strings and genos and alleles:
				genos = [(s,alleles[a]) for s,a in genos]
			yield (c,p,alleles,genos)
		cur2.close()


class TabixReader(pysam.Tabixfile):
	''' A wrapper around pysam.Tabixfile that provides a context and
		allows us to query using 1-based coordinates.
		
		We should request upstream to the Pysam people to implement
		methods __reversed__() and __len__() in pysam.TabixIterator.
		__getitem__ could be a bonus, but prob unnecessary.
	'''
	def __init__(self, inFile, parser=pysam.asTuple()):
		# because of odd Cython weirdness, we don't actually want to call super.__init__ here..
		#super(TabixReader, self).__init__(inFile, parser=parser, mode=mode)
		self.parser = parser
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def close(self):
		if hasattr(super(TabixReader, self), 'close'):
			super(TabixReader, self).close()
		else:
			log.warn("warning: pysam-0.5 lacks a pysam.Tabixfile.close() method.  The input file may not be closed.")
			# NOTE: this does not exist in pysam-0.5
			# I'm not sure if that means it actually doesn't close the file ever.
	def chroms(self):
		return self.contigs
	def get(self, chrom=None, start=None, stop=None, region=None):
		if start!=None:
			start -= 1
		return self.fetch(reference=chrom, start=start, end=stop,
			region=region, parser=self.parser)


if tuple(map(int, pysam.__version__.split('.'))) > (0,5):
	# new versions of pysam return a zero-based position here
	def get_pos_from_vcf_record(vcfrec):
		return vcfrec.pos + 1
else:
	# old versions of pysam return the correct one-based position (same as the VCF file) here
	def get_pos_from_vcf_record(vcfrec):
		return vcfrec.pos


def bytes_to_string(o):
	if type(o) == bytes:
		o = o.decode('utf-8')
	return o

class VcfReader(TabixReader):
	''' Same as TabixReader with a few more perks for VCF files:
		- emit results parsed as pysam VCF rows
		- provide self.chrlens(), a dict mapping chrom name to chrom length
		- provide self.samples(), a list of sample names in order of appearance
		- provide get_range(c,start,stop) and get_snp_genos(c,pos)
	'''
	def __init__(self, inFile, ploidy=1, parser=pysam.asVCF()):
		super(VcfReader, self).__init__(inFile, parser=parser)
		assert ploidy in (1,2)
		self.ploidy = ploidy
		self.clens = []
		self.sample_names = None
		for line in self.header:
			line = bytes_to_string(line)
			if line.startswith('##contig=<ID=') and line.endswith('>'):
				line = line[13:-1]
				c = line.split(',')[0]
				clen = int(line.split('=')[1])
				self.clens.append((c,clen))
			elif line.startswith('#CHROM'):
				row = line.split('\t')
				self.sample_names = row[9:]
		self.clens = dict(self.clens)
		assert self.sample_names
	def samples(self):
		return self.sample_names
	def chrlens(self):
		return self.clens
	def get_positions(self, c=None, start=None, stop=None, region=None):
		for snp in self.get(c,start,stop,region):
			yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), get_pos_from_vcf_record(snp)+len(snp.ref)-1)
	def get_range(self, c=None, start=None, stop=None, region=None, as_strings=True, more=False):
		'''	Read a VCF file (optionally just a piece of it) and return contents
			as an iterator with parsed contents.  Each row is returned as
			a 4-tuple:
				1: chr
				2: pos (int)
				3: list of allele strings (in order)
				4: list of genotypes as 2-tuples:
					haploid: (sample, allele)
					diploid: (sample, [allele, allele])
			If as_strings, the alleles will be actual alleles.  Otherwise,
			alleles will be integers.
			If more is true, a fifth column will be emitted with the pysam VCF object.
		'''
		for snp in self.get(c,start,stop,region):
			alleles = [bytes_to_string(snp.ref)] + bytes_to_string(snp.alt).split(',')
			alleles = [a for a in alleles if a != '.']
			if self.ploidy==1:
				genos = [(self.sample_names[i], int(bytes_to_string(snp[i])[0]))
					for i in range(len(self.sample_names))
					if bytes_to_string(snp[i])[0] != '.']
				if as_strings:
					genos = [(s,alleles[a]) for s,a in genos]
			else:
				genos = [(self.sample_names[i], [int(bytes_to_string(snp[i])[j*2]) for j in range(self.ploidy)])
					for i in range(len(self.sample_names))
					if bytes_to_string(snp[i])[0] != '.']
				if as_strings:
					genos = [(s,[alleles[a] for a in g]) for s,g in genos]
			if more:
				yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), alleles, genos, snp)
			else:
				yield (bytes_to_string(snp.contig), get_pos_from_vcf_record(snp), alleles, genos)
	def get_snp_genos(self, c, p, as_strings=True):
		''' Read a single position from a VCF file and return the genotypes
			as a sample -> allele map.  If there is not exactly one matching
			row in the VCF file at this position (if there are none or multiple)
			then we return an empty map: {}.
		'''
		snps = [x for x in self.get_range(c,p,p,as_strings=as_strings)]
		return len(snps)==1 and dict(snps[0][3]) or {}
	def getFullSequences(self, c, start, stop, samples,
		na='-', refSeq=None, refInvariantsOnly=False, ignoreIndels=False):
		''' chr - chromosome name
			start - start position
			stop - default = start
			samples is a list of samples.  None can be used for the ref sample.
			if refSeq is a string with length = stop-start+1, then use this as the
				base sequence for missing data, otherwise, use the na string.
				if refInvariantsOnly is True, refSeq is only used for invariant
				sites or sites with no entries for any samples.  if False, refSeq
				is used for all missing data.
		'''
		assert 1<=start<=stop
		assert len(na)==1
		
		# get all the VCF records
		vcf_records = [(p-start,alleles,dict(genos)) for chrom,p,alleles,genos
			in self.get_range(c, start, stop, as_strings=True)
			if not ignoreIndels or set(map(len, alleles)) == set([1])]
			
		# Construct a list called "seq" into which we will replace alleles as
		# we discover them.  This is a list, not a string, because each position
		# might be replaced with an arbitrary length string (deletions will be
		# empty strings, insertions will be strings > length 1).  This all gets
		# converted to a string at the end.
		if refSeq:
			# assume refSeq alleles as a baseline for missing data
			assert len(refSeq)==(stop-start+1)
			seq = list(refSeq)
			if refInvariantsOnly:
				# but don't assume refSeq alleles for any known variant sites
				for i,alleles,genos in vcf_records:
					if len(set(genos[s] for s in samples if s in genos))>1:
						for j in range(len(alleles[0])):
							if i+j < len(seq):
								seq[i+j] = na
		else:
			# assume nothing when data is missing
			seq = list(na * (stop-start+1))
		
		for sample in samples:
			assert sample==None or sample in self.samples()
			
			# Make copy of reference sequence
			newseq = [s for s in seq]
	
			# Replace alleles, one site at a time.
			replaceAlleles(sample, newseq, vcf_records)
			
			# Convert list back to string.  This string may or may not be the same
			# length as the original (if ignoreIndels is True, it must be the same
			# length as the original).
			newseq = ''.join(newseq)
			assert len(newseq)==(stop-start+1) or not ignoreIndels
			yield (sample, newseq)


def replaceAlleles(sample,seq,vcf_records):
	''' Replace alleles, one site at a time. '''
	for i,alleles,genos in vcf_records:
		
		# set allele to the DNA sequence we will replace at positions i through i+len(refallele)-1
		if sample==None:
			# caller is asking for the reference sample's sequence
			allele = alleles[0]
			alleles = [allele]
		else:
			# caller is asking for a normal sample
			samp_geno = genos.get(sample)
			if not samp_geno:
				continue
			if isinstance(samp_geno, list):
				log.warn("TO DO: add code to turn hets into IUPAC ambiguity codes (%s %s = %s)." % (i,sample, '/'.join(samp_geno)))
				continue
			allele = samp_geno
		
		# replace portion of sequence with allele
		# NEED BUGFIXES HERE for overlapping VCF records
		if allele == None:
			# nothing here ... is this even possible anymore?
			pass
		elif len(alleles[0])==1:
			# the most common cases: SNP, novar, or indel replacing one ref base (with zero or more bases)
			seq[i] = allele
		else:
			# most complicated case: something replacing multiple ref bases
			# TODO: beware--this is all untested!
			for j in range(max(len(alleles[0]), len(allele))):
				if i+j < len(seq):
					if j<len(allele):
						if j==len(alleles[0])-1:
							# new allele is >= ref length, fill out the rest of the bases
							seq[i+j] = allele[j:]
						else:
							seq[i+j] = allele[j]
					else:
						# new allele is shorter than ref, so delete extra bases
						seq[i+j] = ''
	return seq


class TabixWriter:
	''' A file-writing context that writes tab delimited text to a temp file
		and, upon exit, bgzip compresses to its final destination and tabix
		indexes it.
		tabixPath is only required if
			- header_prefix='' and there are header lines
			- or if you don't have pysam installed
	'''
	def __init__(self, outFile, header_prefix='#',
			col_chr=1, col_start=2, col_stop=2, preset=None, zerobased=False,
			tabixPath=None):
		assert outFile.endswith('.gz')
		assert not os.access(outFile, os.F_OK) or os.access(outFile, os.W_OK)
		self.outFile = outFile
		root_fname = outFile.split('/')[-1].split('.')[0]
		self.tmpFile = util.file.mkstempfname(prefix='tabix-%s-'%root_fname,suffix='.txt')
		self.tmp_f = open(self.tmpFile, 'wt')
		self.params = {'seq_col':col_chr, 'start_col':col_start, 'end_col':col_stop,
			'preset':preset, 'zerobased':zerobased, 'meta_char':header_prefix,
			'line_skip':0}
		self.header = None
		self.n_header = 0
		self.n_rows = 0
		self.tabixPath=tabixPath
		self.last_chr_start = [None,None]
		self.seen_chrs = set()
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close(error=(exc_type!=None))
		return 0
	def close(self, error=False):
		close(self.tmp_f)
		if not error:
			if self.params['meta_char']=='':
				self.params['line_skip'] = self.n_header
			bgzip_index(self.tmpFile, self.outFile, params, tabixPath=self.tabixPath)
			os.unlink(self.tmpFile)
	def write_header(self, header):
		''' header is a list of strings to be joined by tabs '''
		assert not self.n_rows, "cannot write more header again after data has started"
		self.header = header
		self.tmp_f.write(self.params['meta_char']+'\t'.join(header)+'\n')
		self.n_header += 1
	def write_list(self, row):
		''' row is a list of data to be converted to strings and joined by tabs '''
		self.validate(row)
		self.tmp_f.write('\t'.join([x!=None and str(x) or '' for x in row])+'\n')
		self.n_rows += 1
	def write_map(self, row):
		''' row is a dict of data which we key into with the names of the last
			header row. '''
		self.write_list([row.get(h) for h in self.header])
	def validate(self, row):
		c, start, stop = [row[self.params[x]-1] for x in ('seq_col','start_col','end_col')]
		start,stop = (int(start), int(stop))
		assert stop>=start, "error: stop (%d) < start (%s) on chr %s" % (stop,start,c)
		if self.last_chr_start[0]==c:
			assert start>=self.last_chr_start[1], "error: data must be sorted, %s:%d seen after %s:%d" %(c,start,self.last_chr_start[0],self.last_chr_start[1])
		else:
			assert c not in self.seen_chrs, "error: data must be sorted, %s seen in discontinuous order" % c
			self.seen_chrs.add(c)
		self.last_chr_start = [c,start]
		


def bgzip_index(inFile, outFile, params={}, tabixPath=None):
	assert not inFile.endswith('.gz') and outFile.endswith('.gz')
	
	log.debug("compressing with bgzip %s -> %s" % (inFile, outFile))
	if tabixPath:
		cmdline = "%s/bgzip -c %s > %s" % (tabixPath, inFile, outFile)
		assert not os.system(cmdline)
	else:
		pysam.tabix_compress(self.tmpFile, self.outFile, force=True)
	
	log.debug("indexing with tabix: %s" % outFile)
	if tabixPath:
		cmdline = "%s/tabix %s -f" % (tabixPath, outFile)
		if params.get('seq_col')!=None:
			cmdline += ' -s %d' % params['seq_col']
		if params.get('start_col')!=None:
			cmdline += ' -b %d' % params['start_col']
		if params.get('end_col')!=None:
			cmdline += ' -e %d' % params['end_col']
		if params.get('preset')!=None:
			cmdline += ' -p %s' % params['preset']
		if params.get('meta_char')!=None:
			cmdline += ' -c %s' % params['meta_char']
		if params.get('line_skip')!=None:
			cmdline += ' -S %d' % params['line_skip']
		if params.get('zerobased')!=None:
			cmdline += ' -0'
		assert not os.system(cmdline)
	else:
		assert not params.get('line_skip'), "error: pysam does not seem to support this option, even though their documentation talks about it"
		pysam.tabix_index(self.outFile, force=True,
			seq_col=params.get('seq_col'),
			start_col=params.get('start_col'), end_col=params.get('end_col'),
			preset=params.get('preset'), meta_char=params.get('meta_char','#'),
			zerobased=params.get('zerobased',False))
	return outFile
