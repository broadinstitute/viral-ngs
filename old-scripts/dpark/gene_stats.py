#!/usr/bin/env python
'''This script contains a number of tools to pull reports on genes.
Requires python = 2.7, pysam, and BioPython.'''

__author__ = "dpark@broadinstitute.org"
__version__="$Id: gene_stats.py 7870 2014-05-28 14:58:20Z dpark $"[5:-2]
__date__ = "$Date: 2014-05-28 10:58:20 -0400 (Wed, 28 May 2014) $"[7:-2]
__commands__ = []

import argparse, logging
import Bio.Seq, Bio.Alphabet, Bio.SeqIO
import util_cmd, util_files, util_annot, util_vcf

log = logging.getLogger(__name__)
global_tool_paths = {}


def parser_region_intersects():
    parser = argparse.ArgumentParser(
        description='''Add a column to a file of regions that describes the identifiers
of genes that overlap this window.  Take gene list from GFF file.''')
    parser.add_argument("inRegions", help="Input tab file of regions")
    parser.add_argument("inGff", help="Input gene annotation")
    parser.add_argument("outRegions", help="Output tab file with additional column")
    parser.add_argument("--chromInts",
        action="store_true", dest="chromInts",
        default=False,
        help="""Chromosomes are integers in input file and need conversion to strings.  [default: %(default)s]""")
    util_cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',__version__)))
    return parser
def main_region_intersects(args):
    inRegions, inGff, outRegions = (args.inRegions, args.inGff, args.outRegions)

    genedb = util_annot.GeneDb(inGff)
    # read input and write output
    with open(outRegions, 'wt') as outf:
        header = None
        with open(inRegions, 'rt') as inf:
            for line in inf:
                row = line.rstrip('\r\n').split('\t')
                if not header:
                    row.append('genes')
                    header = row
                else:
                    c = row[0]
                    if args.chromInts:
                        c = "Pf3D7_%02d" % int(c)
                    start = int(row[1])
                    stop = int(row[2])
                    genes = [x for x in genedb.getFeaturesInRegion(c,start,stop)]
                    row.append(';'.join(genes))
                outf.write("\t".join(row) + "\n")

    # done
    log.info("done")
    cur.close()
    conn.close()
    return 0
__commands__.append(('region_intersects', main_region_intersects, parser_region_intersects))


def parser_gene_fasta():
    parser = argparse.ArgumentParser(
        description='''Emit a fasta file with gene sequences for every strain based
on the reference sequence substituted with genotype data from a VCF file.
Take gene list from GFF file.''')
    parser.add_argument("inRegions", help="Input tab file of regions")
    parser.add_argument("inGff", help="Input gene annotations")
    parser.add_argument("geneList", help="List of gene IDs")
    parser.add_argument("outFasta", help="Output fasta file of gene sequences")
    parser.add_argument("--sampleList", dest="sampleList",
        help="File with a list of sample names to include (others will be excluded).  One sample per line.",
        default=None)
    parser.add_argument("--type", dest="type", choices=('gene','mRNA','CDS','protein'),
        help="Type of sequence to extract (gene, mRNA, CDS, protein).  [default: %(default)s]",
        default='gene')
    parser.add_argument("--addRefSample",
        action="store_true", dest="addRefSample",
        default=False,
        help="Also output the reference sample's sequence as well.")
    parser.add_argument("--assumeRef", dest="assumeRef",
        help="If specified, this is a fasta file that describes the reference genome.  We will substitute in a reference base any time we have missing data.",
        default=None)
    parser.add_argument("--assumeRefInvariants",
        action="store_true", dest="assumeRefInvariants",
        default=False,
        help="""If false (default), assumeRef will provide a base sequence for
all missing data for each sample.  If true, assumeRef provides a base sequence
for invariant sites only (defined as a site missing from the inVcf or a site
denoted as invariant in the inVcf)""")
    parser.add_argument("--naChar", dest="naChar",
        help="Character to be used for missing data.  Does not apply to type==protein.  [default: %(default)s]",
        default='N')
    parser.add_argument("--ignoreIndels",
        action="store_true", dest="ignoreIndels",
        default=False,
        help="Read only SNPs from the VCF file, ignore any multi-base variations.")
    util_cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',__version__)))
    return parser
def main_gene_fasta(args):
    inVcf, inGff, geneList, outFasta = (args.inVcf, args.inGff, args.geneList, args.outFasta)
    assert len(args.naChar)==1
    type = (args.type=='protein') and 'CDS' or args.type

    # prep inputs
    genedb = util_annot.GeneDb(inGff)
    vcfdb = util_vcf.VcfReader(inVcf)
    with open(geneList, 'rt') as inf:
        geneList = [line.rstrip('\r\n') for line in inf]
        badgenes1 = [gene for gene in geneList if not genedb.getByGeneId(gene, 'gene')]
        if badgenes1:
            log.error("These gene IDs are missing from the GFF file:\n%s" % '\n'.join(badgenes1))
        badgenes2 = [gene for gene in geneList if (gene not in badgenes1) and (not genedb.getByGeneId(gene, type))]
        if badgenes2:
            log.error("These gene IDs are in the GFF file but do not have the %s product requested:\n%s" %  (type, '\n'.join(badgenes2)))
        if badgenes1 or badgenes2:
            return 1
    if args.sampleList:
        with open(args.sampleList, 'rt') as inf:
            sampleList = [line.rstrip('\r\n') for line in inf]
    else:
        sampleList = vcfdb.samples()
    if args.addRefSample:
        sampleList = [None] + sampleList
    refdb = args.assumeRef and Bio.SeqIO.index(open(args.assumeRef, 'rt'), 'fasta') or None
    translator = VerboseTranslator(genedb)

    # start processing and dumping output
    with open(outFasta, 'wt') as outf:
        for gene in geneList:
            log.info("pulling sequence for gene %s across %d individuals" % (gene, len(sampleList)))
            c, start, stop, strand = genedb.getByGeneId(gene)[0]

            # get genetic code for this chromosome if we need to translate
            if args.type=='protein':
                translation_table = translator.get_table(c)

            # get each strains version of this gene (spliced and flipped if needed)
            for sample, seq in mutSpliceSeqs(vcfdb, genedb.getByGeneId(gene, type), sampleList,
                refdb=refdb, naChar=(args.type=='protein') and 'N' or args.naChar,
                refInvariantsOnly=args.assumeRefInvariants, ignoreIndels=args.ignoreIndels):
                if sample==None:
                    sample = 'reference'
                header = "%s:%s:%s" % (gene, sample, type)
                if strand=='-':
                    seq = revcomp(seq)
                if args.type=='protein':
                    seq = translator.translate(seq, translation_table)
                for line in util_files.fastaMaker([(header, seq)]):
                    outf.write(line)

    # done
    genedb.close()
    log.info("done")
    return 0
__commands__.append(('gene_fasta', main_gene_fasta, parser_gene_fasta))


def mutSpliceSeqs(vcf, features, sampleList,
    refdb=None, naChar='N', refInvariantsOnly=True, ignoreIndels=False):
    ''' This function retrieves a specific gene's sequence for a set of samples.
        Each sample will have their own version of the sequence, based on their
        genotypes reported in the vcf file.
        Returns: an iterator of (sample, sequence string) pairs
        Inputs:
            vcf - a VcfReader object
            features - a list of (chrom,start,stop,strand) tuples,
                sorted by start, and all chroms and strands must be identical.
                If len(features)>1, we will splice the resulting sequences together
                at the end on a per-sample basis.
            sampleList - a list of sample names to compute on.  A None entry is
                interpreted as passing the reference sequence unchanged.
            refdb - a Bio.SeqIO.index object called on a reference fasta file.
                Alternatively, this may be None, if no such data is available.
    '''
    samplesToSeqPerExon = [
        dict(vcf.getFullSequences(c,start,stop,sampleList,
                refSeq=(refdb and (refdb[c][start-1:stop].seq) or None),
                refInvariantsOnly=refInvariantsOnly, ignoreIndels=ignoreIndels, na=naChar))
        for c,start,stop,strand in features]
    for sample in sampleList:
        yield (sample, ''.join([samplesToSeq[sample] for samplesToSeq in samplesToSeqPerExon]))

def revcomp(seq):
    return str(Bio.Seq.Seq(seq, Bio.Alphabet.generic_dna).reverse_complement())

class VerboseTranslator:
    def __init__(self, genedb):
        self.genedb = genedb
        self.assumed_chrs = set()
    def get_table(self, c):
        translation_table = self.genedb.getChrInfo(c, 'translation_table')
        if translation_table:
            translation_table = int(translation_table)
        else:
            translation_table = 1
            if c not in self.assumed_chrs:
                log.warn("Assuming standard translation table (standard code 1) for chr %s.  If this is incorrect, you should add a supercontig entry to the GFF file that has a 'translation_table' key-value pair specifying an integer from 1-25 (see NCBI for the list of available genetic codes)." % c)
                self.assumed_chrs.add(c)
        if translation_table != 1:
            log.info("chromosome %s uses non-nuclear translation table %s" % (c,translation_table))
        return translation_table
    def translate(self, dnaseq, translation_table):
        return str(Bio.Seq.Seq(dnaseq, Bio.Alphabet.generic_dna).translate(to_stop=True, table=translation_table))



def translateIfProt(seqIter, features):
    # get strand so we know whether to revcomp or not
    strand = features and features[0][-1] or None

    # revcomp if needed
    for sample, seq in seqIter:
        if strand=='-':
            seq = str(Bio.Seq.Seq(seq, Bio.Alphabet.generic_dna).reverse_complement())
        yield (sample, seq)



def parser_geneReport():
    parser = argparse.ArgumentParser(
        description='''Take a list of gene IDs, a VCF file of SNPs annotated by
snpEff, and emit a tab text file report of all SNPs in those genes.''')
    parser.add_argument("geneList", help="List of gene IDs")
    parser.add_argument("snpEffVcf", help="snpEff-annotated VCF file")
    parser.add_argument("outReport", help="Output tab text report file")
    parser.add_argument("--gff", dest="gff",
        help="If specified, this is a GFF file describing gene coordinates.  This script will execute much faster if this information is provided.",
        default=None)
    parser.add_argument("--diploid",
        action="store_true", dest="diploid",
        default=False,
        help="Report diploid genotypes.  Default is to condense homozygous genotypes into a haploid genotype while passing hets as diploid.")
    util_cmd.common_args(parser, (('tmpDir',None), ('loglevel','INFO'), ('version',__version__)))
    return parser
def main_geneReport(args):
    # read in list of genes
    with open(args.geneList, 'rt') as inf:
        geneList = [x.strip() for x in inf]
    log.info("reporting on %d genes" % len(geneList))

    if args.gff:
        # pull from VCF file based on known coordinates
        with util_annot.GeneDb(args.gff) as genedb:
            badgenes = [gene for gene in geneList if not genedb.getByGeneId(gene)]
            assert not badgenes, "error: genes not present in GFF: %s" % badgenes
            intervals = [genedb.getByGeneId(gene)[0] for gene in geneList]
            intervals = ["%s:%d-%d" % (c,start,stop) for c,start,stop,strand in intervals]
        geneList = None
    else:
        # read through entire VCF file and filter manually based on gene names (slow)
        intervals = [None]
        geneList = set(geneList)

    # execute
    vcf_to_text_report(args.snpEffVcf, intervals, args.outReport, args.diploid and 2 or 1, gene_list=geneList)

    # done
    log.info("done")
    return 0
__commands__.append(('gene_report', main_geneReport, parser_geneReport))


def parser_regionReport():
    parser = argparse.ArgumentParser(
        description='''Take a VCF file (optionally annotated by snpEff),
specified regions, and emit a tab text file report of all SNPs in these windows.
If the interval list is simply 'all', we will emit the entire VCF.''')
    parser.add_argument("inVcf", help="Input VCF file of variants")
    parser.add_argument("intervals", help="Intervals to report on: chr:start-stop,chr:start-stop,...")
    parser.add_argument("outReport", help="Output tab text report file")
    parser.add_argument("--diploid",
        action="store_true", dest="diploid",
        default=False,
        help="Report diploid genotypes.  Default is to condense homozygous genotypes into a haploid genotype while passing hets as diploid.")
    util_cmd.common_args(parser, (('tmpDir',None), ('loglevel','INFO'), ('version',__version__)))
    return parser
def main_regionReport(args):
    if args.intervals=='all':
        intervals = [None]
    else:
        intervals = args.intervals.split(',')
    vcf_to_text_report(args.inVcf, intervals, args.outReport, args.diploid and 2 or 1)
    log.info("done")
    return 0
__commands__.append(('region_report', main_regionReport, parser_regionReport))


def vcf_to_text_report(inVcf, intervals, outReport, ploidy=1, gene_list=None, all_genos=True):
    ''' Read a VCF file at the specified intervals and emit each row in a somewhat
        modified tab text format for easy human readability.
        If intervals is [None], we will emit the whole VCF.
        If gene_list is specified (a set or dict), we will filter out any row
        that does not specify one of the genes in gene_list in a snpEff-formatted
        fashion in the INFO field of the VCF.
        If all_genos is false, we will omit the per-sample genotypes from the output.
    '''
    assert ploidy==1, "non-haploid not working yet"

    # write output file
    with open(outReport, 'wt') as outf:
        pre_header = ['chr','pos','gene','desc','effect','impact',
            'protein_pos','residue_ref','residue_alt',
            'allele_minor','allele_ref','allele_alt','call_rate','mac','maf']

        with util_vcf.VcfReader(inVcf) as vcf:  #, ploidy=ploidy -- why doesn't this work?
            # get header
            indivs = vcf.samples()
            log.info("reporting on %d individuals" % len(indivs))
            header = all_genos and (pre_header+indivs) or pre_header
            outf.write('#' + '\t'.join(header) + '\n')

            # grab each interval
            for interval in intervals:
                log.info("fetching interval %s" % interval)
                # emit per snp
                for c,p,alleles,genos,snp in vcf.get_range(region=interval, more=True):
                    if len(alleles)>1:
                        # prep all the stats/metadata on this SNP
                        eff = util_annot.parse_eff(c,p,snp.info,required=False)
                        stats = {
                            'chr':c,
                            'pos':p,
                            'allele_ref':alleles[0],
                            'allele_alt':','.join(alleles[1:]),
                            'effect':eff[0],
                            'impact':eff[1],
                            'gene':eff[2],
                            'desc':eff[3],
                            'protein_pos':eff[4],
                            'residue_ref':eff[5],
                            'residue_alt':eff[6],
                        }
                        if gene_list and stats['gene'] not in gene_list:
                            continue

                        # parse genos
                        genos = dict(genos)
                        if ploidy==2:
                            for k in genos.keys():
                                genos[k] = '/'.join(genos[k])

                        # calc call rate and maf
                        stats['call_rate'] = float(len(genos)) / len(indivs)
                        counts = util_vcf.calc_maf(genos.values(),ploidy=ploidy)
                        stats['allele_minor'] = counts['a_minor']
                        stats['mac'] = counts['mac']
                        stats['maf'] = counts['maf']

                        # output to table
                        if all_genos:
                            outf.write('\t'.join(map(str,
                                [stats[x] for x in pre_header] + [genos.get(s,'') for s in indivs])) + '\n')
                        else:
                            outf.write('\t'.join(map(str,
                                [stats[x] for x in pre_header])) + '\n')
    return


if __name__ == '__main__':
    util_cmd.main_argparse(__commands__, __version__, global_tool_paths, __doc__)

