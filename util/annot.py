'''A few methods for dealing with gene annotation from snpEff.

Requires python >= 2.6'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import sqlite3, itertools, urllib, logging, re, os
import util.files

log = logging.getLogger(__name__)

def unique(items):
	''' Return unique items in the same order as seen in the input. '''
	seen = set()
	for i in items:
		if i not in seen:
			seen.add(i)
			yield i

class SnpAnnotater:
	''' Add annotations to snps based on a snpEff-annotated VCF file.
	'''
	def __init__(self, snpEffVcf=None, snpIterator=None):
		self.snpIterator = snpIterator
		self.dbFile = util_files.mkstempfname(prefix='SnpAnnotater-', suffix='.db')
		self.conn = sqlite3.connect(self.dbFile, isolation_level='DEFERRED')
		self.cur = self.conn.cursor()
		self.cur.execute("""create table annot (
			chr not null,
			pos integer not null,
			allele_ref not null,
			allele_alt not null,
			effect not null,
			impact not null,
			gene_id,
			gene_name,
			protein_pos integer,
			residue_ref,
			residue_alt
		)""")
		self.cur.execute("create index idx_annot on annot(chr,pos)")
		if snpEffVcf:
			self.loadVcf(snpEffVcf)
	def loadVcf(self, snpEffVcf):
		#log.info("reading in snpEff VCF file: %s" % snpEffVcf)
		with util_files.open_or_gzopen(snpEffVcf, 'rt') as inf:
			ffp = util_files.FlatFileParser(inf)
			try:
				imap = hasattr(itertools, 'imap') and itertools.imap or map  #py2 & py3 compatibility
				ifilter = hasattr(itertools, 'ifilter') and itertools.ifilter or filter  #py2 & py3 compatibility
				self.cur.executemany("""insert into annot (chr,pos,allele_ref,allele_alt,
					effect,impact,gene_id,gene_name,protein_pos,residue_ref,residue_alt)
					values (?,?,?,?,?,?,?,?,?,?,?)""",
					imap(lambda row:
						[row['CHROM'], int(row['POS']), row['REF'], row['ALT']]
						+ parse_eff(row['CHROM'], row['POS'], row['INFO']),
						ifilter(lambda r: r['ALT'] != '.', ffp)))
			except Exception as e:
				log.exception("exception processing file %s line %s" % (snpEffVcf, ffp.line_num))
				raise
			self.cur.execute("select chr,pos from annot group by chr,pos having count(*)>1")
			dupes = [(c,p) for c,p in self.cur]
			if dupes:
				log.info("deleting annotation for %d duplicate positions: %s" % (len(dupes), ', '.join(['%s:%s'%(c,p) for c,p in dupes])))
				self.cur.executemany("delete from annot where chr=? and pos=?", dupes)
			self.conn.commit()
	def __iter__(self):
		assert self.snpIterator
		for snp in self.snpIterator:
			yield self.annotate(snp)
	def annotate(self, row):
		self.cur.execute("""select effect,impact,gene_id,gene_name,protein_pos,
			allele_ref,allele_alt,residue_ref,residue_alt
			from annot where chr=? and pos=?""", [row['chr'], int(row['pos'])])
		x = self.cur.fetchone()
		if x != None:
			row['effect'],row['impact'],row['gene_id'],row['gene_name'],row['protein_pos'],row['allele_ref'],row['allele_alt'],row['residue_ref'],row['residue_alt'] = x
			row['alleles'] = '/'.join((row['allele_ref'],row['allele_alt']))
			if row['residue_alt']:
				row['residues'] = '/'.join((row['residue_ref'],row['residue_alt']))
			else:
				row['residues'] = row['residue_ref']
		else:
			row['effect'] = 'UNKNOWN'
			row['impact'] = 'UNKNOWN'
		return row
	def new_fields(self):
		return ('effect', 'impact', 'gene_id', 'gene_name',	'protein_pos', 'alleles', 'residues')
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def close(self):
		self.cur.close()
		self.conn.close()
		os.unlink(self.dbFile)



def parse_eff(chrom, pos, info, required=True):
	try:
		impact_rank = {'HIGH':0,'MODERATE':1,'LOW':2,'MODIFIER':3}
		infos = [x for x in info.split(';') if x.startswith('EFF=')]
		assert len(infos)<=1
		if not infos:
			assert not required
			return ['', '', '', '', '', '', '']
		info = infos[0]
		effs = info[4:].split(',')
		out = []
		for eff in effs:
			assert eff.endswith(')')
			effect, other = eff[:-1].split('(')
			other = other.split('|')
			assert len(other)>=10
			impact = other[0]
			gene_id = other[8]
			assert not gene_id or (gene_id.endswith('-1') and gene_id.startswith('rna_'))
			if gene_id:
				gene_id = gene_id[4:-2]
			if gene_id=='PF14_0620' or gene_id=='PF3D7_1465300':
				gene_name = 'tRNA 3-trailer sequence RNase, putative'
			else:
				try:
					gene_name = urllib.unquote_plus(other[5]).encode('ascii')
				except UnicodeDecodeError as e:
					log.error("error at %s:%s decoding the string '%s'" % (chrom, pos, other[5]))
					raise
			aa_chg = other[3]
			if aa_chg:
				if effect.startswith('SYNON'):
					aas = (aa_chg[0], '')
					codon = int(aa_chg[1:])
				elif effect.startswith('NON_SYNON') or effect.startswith('START_') or effect.startswith('STOP_') or effect.startswith('CODON_'):
					mo = re.match(r'^(\D+)(\d+)(\D+)$', aa_chg)
					assert mo, "unrecognized coding change format for %s (%s)" % (effect, aa_chg)
					aas = (mo.group(1), mo.group(3))
					codon = int(mo.group(2))
				elif effect=='FRAME_SHIFT':
					mo = re.match(r'^(\D*)(\d+)(\D*)$', aa_chg)
					assert mo, "unrecognized coding change format for %s (%s)" % (effect, aa_chg)
					aas = ('','')
					codon = int(mo.group(2))
				else:
					assert 0, "unrecognized effect type (%s) for variant with coding change (%s)" % (effect, aa_chg)
			else:
				aas = ('','')
				codon = ''
			out.append([impact_rank[impact], effect, impact, gene_id, gene_name,
				codon, aas[0], aas[1]])
	
		if len(out)>1:
			out.sort()
			if out[0][0] == out[1][0]:
				#log.debug("SNP found with multiple effects of the same impact level: %s:%s - %s" % (chrom, pos, info))
				#assert out[0][2] in ('MODIFIER','LOW'), "Error: SNP found with multiple effects of the same impact level"
				out = [[';'.join(unique([str(o[i]) for o in out])) for i in range(len(out[0]))]]
		eff = out[0][1:]
		return eff
		
	except Exception as e:
		log.exception("exception parsing snpEff on row %s:%s - %s" % (chrom, pos, info))
		raise
	except Error as e:
		log.error("error parsing snpEff on row %s:%s - %s" % (chrom, pos, info))
		raise



class DbConnection:
	def __init__(self, dbFile=None):
		if dbFile==None:
			dbFile = util_files.mkstempfname(suffix='.db')
		self.conn = sqlite3.connect(dbFile, isolation_level='DEFERRED')
		assert self.conn.isolation_level
		self.conn.text_factory = sqlite3.OptimizedUnicode
		self.cur = self.conn.cursor()
		self.cur.execute("PRAGMA foreign_keys=ON")
		self.cur.execute("PRAGMA foreign_keys")
		fk = self.cur.fetchone()
		log.debug("SQLite version: %s" % sqlite3.sqlite_version)
		log.debug("SQLite foreign key support: %s" % ((fk and fk[0]) and 'true' or 'false'))
		self.start()
	def start(self):
		pass
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()
		return 0
	def close(self):
		if not self.conn:
			return
		self.conn.commit()
		self.cur.close()
		self.conn.close()
		self.conn = None
		self.cur = None

class GeneDb(DbConnection):
	def __init__(self, gffFile=None):
		DbConnection.__init__(self)
		if gffFile:
			self.loadGff(gffFile)
	def start(self):
		DbConnection.start(self)
		self.cur.execute("""create table if not exists gene (
			id text primary key,
			chr text not null,
			start integer not null,
			stop integer not null,
			strand text not null,
			name text)""")
		self.cur.execute("create index idx_gene_coords on gene(chr,start,stop)")
		self.cur.execute("""create table if not exists child (
			id text not null,
			gene_id text not null,
			start integer not null,
			stop integer not null,
			type text not null,
			foreign key (gene_id) references gene(id))""")
		self.cur.execute("create index idx_child_gene on child(gene_id,type,start)")
		self.cur.execute("create index idx_child_idx on child(id,type)")
		self.cur.execute("""create table if not exists chr_info (
			chr text not null,
			key text not null,
			val text not null)""")
		self.cur.execute("create unique index idx_chr_info on chr_info(chr,key)")
	def loadGff(self, inGff):
		log.info("loading GFF %s" % inGff)
		with util_files.open_or_gzopen(inGff, 'rt') as inf:
			line_num = 0
			for line in inf:
				line_num += 1
				row = line.rstrip('\r\n').split('\t')
				if line.startswith('#'):
					continue
				c = row[0]
				type = row[2]
				start = int(row[3])
				stop = int(row[4])
				strand = row[6]
				info = row[8].split(';')
				assert info[0].startswith('ID=')
				id = info[0][3:]
				if type == 'supercontig':
					info = dict([x.split('=') for x in info])
					self.cur.executemany("insert into chr_info (chr,key,val) values (?,?,?)",
						[(c,k,v) for k,v in info.items()])
				elif type == 'gene':
					assert info[1].startswith('Name=')
					name = info[1][5:]
					if id=='PF14_0620' or id=='PF3D7_1465300':
						name = 'tRNA 3-trailer sequence RNase, putative'
					else:
						name = urllib.unquote_plus(name).encode('ascii')
					self.cur.execute("insert into gene (chr,start,stop,strand,id,name) values (?,?,?,?,?,?)",
						[c,start,stop,strand,id,name])
				elif type == 'mRNA':
					gene_id = [x[7:] for x in info if x.startswith('Parent=')]
					assert len(gene_id)==1
					gene_id = gene_id[0]
					self.cur.execute("select strand from gene where id=?", [gene_id])
					s = self.cur.fetchone()
					assert s and s[0]==strand
					self.cur.execute("insert into child (id,gene_id,start,stop,type) values (?,?,?,?,?)",
						[id,gene_id,start,stop,type])
				elif type in ('CDS','exon'):
					mrna_id = [x[7:] for x in info if x.startswith('Parent=')]
					assert len(mrna_id)==1
					self.cur.execute("select gene_id from child where id=?", mrna_id)
					gene_id = self.cur.fetchone()
					if gene_id:
						#assert gene_id, "error on line %d: never heard of mRNA %s when loading %s %s" % (line_num, mrna_id[0], type, id)
						gene_id = gene_id[0]
						self.cur.execute("insert into child (id,gene_id,start,stop,type) values (?,?,?,?,?)",
							[id,gene_id,start,stop,type])
		self.conn.commit()
		counts={}
		self.cur.execute("select count(*) from gene")
		counts['gene'] = self.cur.fetchone()[0]
		self.cur.execute("select count(*) from child where type='mRNA'")
		counts['mRNA'] = self.cur.fetchone()[0]
		self.cur.execute("select count(*) from child where type='exon'")
		counts['exon'] = self.cur.fetchone()[0]
		self.cur.execute("select count(*) from child where type='CDS'")
		counts['CDS'] = self.cur.fetchone()[0]
		self.cur.execute("select count(distinct g.id) from gene g join child m on g.id=m.gene_id and m.type='mRNA'")
		counts['gene_with_mRNA'] = self.cur.fetchone()[0]
		self.cur.execute("select count(distinct g.id) from gene g join child c on g.id=c.gene_id and c.type='CDS'")
		counts['gene_with_CDS'] = self.cur.fetchone()[0]
		log.info("GFF file loaded with %d genes, %d mRNAs, %d exons, and %d CDS regions" %
			(counts['gene'], counts['mRNA'], counts['exon'], counts['CDS']))
		log.info("%d genes have mRNAs, %d genes have CDS regions" %
			(counts['gene_with_mRNA'], counts['gene_with_CDS']))
	def getChrInfo(self, c, k=None):
		if k==None:
			self.cur.execute("select key,val from chr_info where chr=?", [c])
			return dict([(k,v) for k,v in self.cur])
		else:
			self.cur.execute("select val from chr_info where chr=? and key=?", [c,k])
			row = self.cur.fetchone()
			if row:
				return row[0]
			else:
				return None
	def getByGeneId(self, gene_id, type='gene'):
		self.cur.execute("select chr,start,stop,strand from gene where id=?", [gene_id])
		feature = self.cur.fetchone()
		if feature==None:
			out = []
		elif type!='gene':
			chr,start,stop,strand = feature
			self.cur.execute("select start,stop from child where gene_id=? and type=? order by start", [gene_id,type])
			out = [(chr,x[0],x[1],strand) for x in self.cur]
		else:
			out = [feature]
		return out
	def getFeaturesInRegion(self, c, start, stop, type='gene'):
		if type=='gene':
			self.cur.execute("select id from gene where chr=? and stop>=? and start<=? order by start", [c,start,stop])
		else:
			self.cur.execute("select id from child where chr=? and stop>=? and start<=? and type=? order by start", [c,start,stop,type])
		for x in self.cur:
			yield x[0]
	def getGeneNameById(self, gene_id):
		self.cur.execute("select name from gene where id=?", [gene_id])
		x = self.cur.fetchone()
		return x!=None and x[0] or None

