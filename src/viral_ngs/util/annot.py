'''A few methods for dealing with gene annotation from snpEff.'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import sqlite3, itertools, urllib, logging, re, os
import util.file, util.misc

log = logging.getLogger(__name__)

class SnpAnnotater:
    ''' Add annotations to snps based on a snpEff-annotated VCF file.
    '''
    def __init__(self, snpEffVcf=None, snpIterator=None):
        self.snpIterator = snpIterator
        self.dbFile = util.file.mkstempfname(prefix='SnpAnnotater-', suffix='.db')
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
        with util.file.open_or_gzopen(snpEffVcf, 'rt') as inf:
            ffp = util.file.FlatFileParser(inf)
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
        return ('effect', 'impact', 'gene_id', 'gene_name', 'protein_pos', 'alleles', 'residues')
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
                out = [[';'.join(util.misc.unique([str(o[i]) for o in out])) for i in range(len(out[0]))]]
        eff = out[0][1:]
        return eff

    except Exception as e:
        log.exception("exception parsing snpEff on row %s:%s - %s" % (chrom, pos, info))
        raise
    except Error as e:
        log.error("error parsing snpEff on row %s:%s - %s" % (chrom, pos, info))
        raise


