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


class TabixReader(pysam.Tabixfile):
    ''' A wrapper around pysam.Tabixfile that provides a context and
        allows us to query using 1-based coordinates.

        We should request upstream to the Pysam people to implement
        methods __reversed__() and __len__() in pysam.TabixIterator.
        __getitem__ could be a bonus, but prob unnecessary.
    '''
    def __init__(self, inFile, parser=pysam.asTuple()):
        # because of odd Cython weirdness, we don't actually want to call super.__init__ here..
        #super(TabixReader, self).__init__(inFile, parser=parser)
        self.parser = parser
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return 0
    def close(self):
        super(TabixReader, self).close()
    def chroms(self):
        return self.contigs
    def get(self, chrom=None, start=None, stop=None, region=None):
        if start!=None:
            start -= 1
        return self.fetch(reference=chrom, start=start, end=stop,
            region=region, parser=self.parser)


def get_pos_from_vcf_record(vcfrec):
    # new versions of pysam return a zero-based position here
    return vcfrec.pos + 1

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
        ''' Read a VCF file (optionally just a piece of it) and return contents
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

