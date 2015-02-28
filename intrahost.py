#!/usr/bin/env python
'''This script contains a number of utilities for intrahost variant calling
and annotation for viral genomes.
'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org, swohl@broadinstitute.org, irwin@broadinstitute.org"
__commands__ = []

import argparse, logging, itertools, re, os
import Bio.AlignIO, Bio.SeqIO, Bio.Data.IUPACData
import pysam
import util.cmd, util.file, util.vcf, util.misc
from util.misc import mean, median
from interhost import CoordMapper
from tools.vphaser2 import Vphaser2Tool

log = logging.getLogger(__name__)


#  ============= class AlleleFieldParser =================

class AlleleFieldParser(object) :
    """
    Class for parsing fields in the allele columns of vphaser_one_sample output
    (corresponding to the SNP_or_LP_Profile columns in the V-Phaser 2 output).
    """
    def __init__(self, field) :
        """Input is the string stored in one of the allele columns."""
        words = field.split(':')
        self.allele = words[0]
        self.libCounts = [[int(words[1]), int(words[2])]]
        self.libBiasPval = 1.0
    
    def get_strand_counts(self) :
        "Return [# forward reads (all libs), # reverse reads (all libs)]"
        return map(sum, zip(*self.libCounts))
    
    def get_allele(self) :
        """ Return allele:
                A, C, G, or T for SNVs,
                Dn with n > 0 for deletions,
                Ibases where bases is a string of two or more bases for inserts
        """
        return self.allele
    
    def get_lib_counts(self) :
        """Yield [# forward reads, # reverse reads] for each library, in order
            of read groups in BAM file.
        """
        for counts in self.libCounts :
            yield counts

    def get_lib_bias_pval(self) :
        "Return a p-value on whether there is a library bias for this allele."
        return self.libBiasPval


#  ========== vphaser_one_sample =================

defaultMinReads = 5
defaultMaxBias = 10

def vphaser_one_sample(inBam, outTab, vphaserNumThreads = None,
                       minReadsEach = None, maxBias = None) :
    ''' Input: a single BAM file, representing reads from one sample, mapped to
            its own consensus assembly. It may contain multiple read groups and 
            libraries.
        Output: a tab-separated file with no header containing filtered
            V Phaser-2 output variants with additional columns:
                sequence/chrom name, # libraries, chi-sq for library discordance
    '''
    if minReadsEach != None :
        assert minReadsEach > 0, 'minReadsEach must be at least 1.'
    variantIter = Vphaser2Tool().iterate(inBam, vphaserNumThreads)
    filteredIter = filter_strand_bias(variantIter, minReadsEach, maxBias)
    libraryFilteredIter = filter_library_bias(filteredIter)
    with open(outTab, 'wt') as outf :
        for row in libraryFilteredIter :
            outf.write('\t'.join(row) + '\n')

def filter_strand_bias(isnvs, minReadsEach = None, maxBias = None) :
    ''' Take an iterator of V-Phaser output (plus chromosome name prepended)
        and perform hard filtering for strand bias
    '''
    if minReadsEach == None :
        minReadsEach = defaultMinReads
    if maxBias == None :
        maxBias = defaultMaxBias
    for row in isnvs:
        front = row[:7]
        acounts = [x.split(':') for x in row[7:]]
        acounts = list([(a,f,r) for a,f,r in acounts
            if int(f)>=minReadsEach and int(r)>=minReadsEach
            and maxBias >= (float(f)/float(r)) >= 1.0/maxBias])
        if len(acounts) > 1:
            acounts = list(reversed(sorted((int(f)+int(r),a,f,r) for a,f,r in acounts)))
            mac = sum(n for n,a,f,r in acounts[1:])
            tot = sum(n for n,a,f,r in acounts)
            back = [':'.join([a,f,r]) for n,a,f,r in acounts]
            front[2] = acounts[1][1]
            front[3] = acounts[0][1]
            front[6] = '%.6g' % (100.0*mac/tot)
            yield front + back

def filter_library_bias(isnvs) :
    ''' Filter variants based on library bias. For ones that pass the filter
            add fields with the number of libraries and a bias p-value.
        NOT YET IMPLEMENTED!
    '''
    for row in isnvs :
        strNlibs = ''   # To be filled in in future
        strLibBias = '' # To be filled in in future
        row = row[:7] + [strNlibs, strLibBias] + row[7:]
        yield row

def parser_vphaser_one_sample(parser = argparse.ArgumentParser()) :
    parser.add_argument("inBam",
        help = "Input Bam file representing reads from one sample, mapped to "
               "its own consensus assembly. It may contain multiple read "
               "groups and libraries.")
    parser.add_argument("outTab", help = "tab-separated headerless output file.")
    parser.add_argument("--vphaserNumThreads", type = int, default = None,
        help="Number of threads in call to V-Phaser 2.")
    parser.add_argument("--minReadsEach", type = int, default = None,
        help = "Minimum number of reads on each strand (default: %s). Must be "
               "at least 1." %
        defaultMinReads)
    parser.add_argument("--maxBias", type = int, default = None,
        help = "Maximum allowable ratio of number of reads on the two strands "
               "(default: %s)." % defaultMaxBias)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, vphaser_one_sample, split_args = True)
    return parser
__commands__.append(('vphaser_one_sample', parser_vphaser_one_sample))


#  ========== tabfile_values_rename =================

def tabfile_values_rename(inFile, mapFile, outFile, col=0):
    ''' Take input tab file and copy to an output file while changing
        the values in a specific column based on a mapping file.
        The first line will pass through untouched (it is assumed to be
        a header).
    '''
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
def parser_tabfile_rename(parser=argparse.ArgumentParser()):
    parser.add_argument("inFile", help="Input flat file")
    parser.add_argument("mapFile",
        help="""Map file.  Two-column headerless file that maps input values to
        output values.  This script will error if there are values in inFile that do
        not exist in mapFile.""")
    parser.add_argument("outFile", help="Output flat file")
    parser.add_argument("--col_idx", dest="col", type=int,
        help="""Which column number to replace (0-based index). [default: %(default)s]""",
        default=0)
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, tabfile_values_rename, split_args=True)
    return parser
__commands__.append(('tabfile_rename', parser_tabfile_rename))


#  ========== merge_to_vcf ===========================

def merge_to_vcf(refFasta, outVcf, samples, isnvs, assemblies):
    ''' Convert vPhaser2 parsed filtered output text file into VCF format.
    '''
    
    # setup
    if not (len(samples) == len(isnvs) == len(assemblies)):
        raise LookupError("samples, isnvs, and assemblies must have the same number of elements")
    samp_to_fasta = dict(zip(samples, assemblies))
    samp_to_isnv = dict(zip(samples, isnvs))
    samp_to_cmap = dict((s, CoordMapper(genome, refFasta))
        for s, genome in samp_to_fasta.items())
    with open(refFasta, 'rU') as inf:
        ref_chrlens = list((seq.id, len(seq)) for seq in Bio.SeqIO.parse(inf, 'fasta'))
    
    if outVcf.endswith('.vcf.gz'):
        tmpVcf = util.file.mkstempfname('.vcf')
    elif outVcf.endswith('.vcf'):
        tmpVcf = outVcf
    else:
        raise ValueError("outVcf must end in .vcf or .vcf.gz")
    
    with open(tmpVcf, 'wt') as outf:
    
        # write header
        outf.write('##fileformat=VCFv4.1\n')
        outf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outf.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        for c, clen in ref_chrlens:
            outf.write('##contig=<ID=%s,length=%d>\n' % (c, clen))
        outf.write('##reference=file://%s\n' % refFasta)
        header = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples
        outf.write('#'+'\t'.join(header)+'\n')

        # one reference chrom at a time
        with open(refFasta, 'rU') as inf:
            for ref_sequence in Bio.SeqIO.parse(inf, 'fasta'):

                # read in all iSNVs for this chrom and map to reference coords
                data = []
                for s in samples:
                    for row in util.file.read_tabfile(samp_to_isnv[s]):
                        s_chrom = samp_to_cmap[s].mapBtoA(ref_sequence.id, None)[0]
                        if row[0] == s_chrom:
                            row = {'sample':s, 'CHROM':ref_sequence.id,
                                's_chrom':s_chrom, 's_pos':int(row[1]),
                                's_alt':row[2], 's_ref': row[3],
                                'n_libs':row[7], 'lib_bias': row[8],
                                'alleles':list(x.split(':') for x in row[9:] if x)}
                            # make a sorted allele list
                            row['allele_counts'] = list(sorted(
                                [(a,int(f)+int(r)) for a,f,r in row['alleles']],
                                key=(lambda x:x[1]), reverse=True))
                            # reposition vphaser deletions minus one to be consistent with
                            # VCF conventions
                            if row['s_alt'].startswith('D'):
                                for a,n in row['allele_counts']:
                                    assert a[0] in ('D','i')
                                row['s_pos'] = row['s_pos']-1
                            # map position back to reference coordinates
                            row['ref_range'] = samp_to_cmap[s].mapAtoB(s_chrom, row['s_pos'])[1]
                            row['POS'] = row['ref_range'] if (type(row['ref_range'])==int) else row['ref_range'][0]
                            data.append(row)
            
                # sort all iSNVs (across all samples) and group by position
                data = sorted(data, key=(lambda row: row['POS']))
                data = itertools.groupby(data, lambda row: row['POS'])
            
                # process one reference position at a time
                for pos, rows in data:
                    rows = list(rows)
                
                    # define the length of this variation based on the largest deletion
                    end = pos
                    for row in rows:
                        for a,n in row['allele_counts']:
                            if a.startswith('D'):
                                # end of deletion in sample's coord space
                                local_end = row['s_pos']+int(a[1:])
                                # end of deletion in reference coord space
                                ref_end = samp_to_cmap[row['sample']].mapAtoB(row['s_chrom'], local_end)
                                if not type(ref_end) == int:
                                    ref_end = ref_end[1]
                                end = max(end, ref_end)
                
                    # find reference allele and consensus alleles
                    refAllele = str(ref_sequence[pos-1:end].seq)
                    consAlleles = {}
                    for s in samples:
                        cons = Bio.SeqIO.index(samp_to_fasta[s], 'fasta')[samp_to_cmap[s].mapBtoA(ref_sequence.id, None)[0]]
                        cons_start = samp_to_cmap[s].mapBtoA(ref_sequence.id, pos)[1]
                        cons_start = cons_start if type(cons_start)==int else cons_start[0]
                        cons_stop  = samp_to_cmap[s].mapBtoA(ref_sequence.id, end)[1]
                        cons_stop  = cons_stop if type(cons_stop)==int else cons_stop[1]
                        allele = str(cons[cons_start-1:cons_stop].seq).upper()
                        if all(a in set(('A','C','T','G')) for a in allele):
                            consAlleles[s] = allele
                        else:
                            log.warn("dropping unclean consensus for %s at %s-%s: %s" % (s, pos, end, allele))
                    
                    # define genotypes and fractions
                    iSNVs = {}
                    for s in samples:
                        
                        # get all rows for this sample and merge allele counts together
                        acounts = dict(itertools.chain.from_iterable(row['allele_counts']
                            for row in rows if row['sample'] == s))
                        if 'i' in acounts and 'd' in acounts:
                            # this sample has both insertions and deletions at the same spot!
                            # average the two (discordant) reference allele counts and
                            # drop one of them (so we have only one reference allele
                            acounts['i'] = int(round((acounts['i'] + acounts['d'])/2.0,0))
                            del acounts['d']
                        
                        if acounts:
                            # we have iSNV data on this sample
                            consAllele = consAlleles[s]
                            tot_n = sum(acounts.values())
                            iSNVs[s] = {}
                            for a,n in acounts.items():
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
                                        log.warn("vPhaser and assembly pipelines mismatch at %s:%d %s - consensus %s, vPhaser %s" % (ref_sequence.id, pos, s, consAllele[0], a))
                                    a = a + consAllele[1:]
                                assert a and a==a.upper()
                                iSNVs[s][a] = f
                            if all(len(a)==1 for a in iSNVs[s].keys()):
                                assert consAllele in iSNVs[s].keys()
                        elif s in consAlleles:
                            # there is no iSNV data for this sample, so substitute the consensus allele
                            iSNVs[s] = {consAlleles[s]:1.0}

                    # get unique alleles list for this position
                    # allele list should start with the reference allele
                    # and should contain all alternate alleles sorted in
                    # descending order, first by consensus-level allele
                    # frequency, second by intrahost read frequency summed
                    # over the population, third by the allele string itself
                    alleles_cons = [a for a,n in sorted(util.misc.histogram(consAlleles.values()).items(), key=lambda x:x[1], reverse=True) if a!=refAllele]
                    alleles_isnv = list(itertools.chain.from_iterable([iSNVs[s].items() for s in samples if s in iSNVs]))
                    alleles_isnv2 = []
                    for a in set(a for a,n in alleles_isnv):
                        counts = list(x[1] for x in alleles_isnv if x[0]==a)
                        alleles_isnv2.append((len(counts),sum(counts),a))
                    alleles_isnv = list(allele for n_samples, n_reads, allele in reversed(sorted(alleles_isnv2)))
                    alleles = list(util.misc.unique([refAllele] + alleles_cons + alleles_isnv))
                    
                    # map alleles from strings to numeric indexes
                    assert len(alleles)>1
                    alleleMap = dict((a,i) for i,a in enumerate(alleles))
                    genos = [str(alleleMap.get(consAlleles.get(s),'.')) for s in samples]
                    freqs = [(s in iSNVs) and ','.join(map(str, [iSNVs[s].get(a,0.0) for a in alleles[1:]])) or '.' for s in samples]

                    # prepare output row and write to file
                    out = [ref_sequence.id, pos, '.', alleles[0], ','.join(alleles[1:]), '.', '.', '.', 'GT:AF']
                    out = out + list(map(':'.join, zip(genos, freqs)))
                    outf.write('\t'.join(map(str, out))+'\n')
    
    # compress output if requested
    if outVcf.endswith('.vcf.gz'):
        pysam.tabix_compress(tmpVcf, outVcf, force=True)
        pysam.tabix_index(outVcf, force=True, preset='vcf')
        os.unlink(tmpVcf)

def parser_merge_to_vcf(parser=argparse.ArgumentParser()):
    parser.add_argument("refFasta",
        help="""The target reference genome. outVcf will use these
            chromosome names, coordinate spaces, and reference alleles""")
    parser.add_argument("outVcf",
        help="Output VCF file containing all variants")
    parser.add_argument("--samples", nargs='+', required=True,
        help="A list of sample names")
    parser.add_argument("--isnvs", nargs='+', required=True,
        help="""A list of file names from the output of vphaser_one_sample
            These must be in the SAME ORDER as samples.""")
    parser.add_argument("--assemblies", nargs='+', required=True,
        help="""a list of consensus fasta files that were used as the
            per-sample reference genomes for generating isnvs.
            These must be in the SAME ORDER as samples.""")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, merge_to_vcf, split_args=True)
    return parser
__commands__.append(('merge_to_vcf', parser_merge_to_vcf))


#  ===================================================

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
    '''Compute the Fws statistic on iSNV data. See Manske, 2012 (Nature)'''
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

def parser_Fws(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outVcf", help="Output VCF file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, add_Fws_vcf, split_args=True)
    return parser
__commands__.append(('Fws', parser_Fws))


#  ===================================================

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

def parser_iSNV_table(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outFile", help="Output text file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_iSNV_table)
    return parser
def main_iSNV_table(args):
    '''Convert VCF iSNV data to tabular text'''
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
__commands__.append(('iSNV_table', parser_iSNV_table))


#  ===================================================

def iSNP_per_patient(table, agg_fun=median):
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
def parser_iSNP_per_patient(parser=argparse.ArgumentParser()):
    parser.add_argument("inFile", help="Input text file")
    parser.add_argument("outFile", help="Output text file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_iSNP_per_patient)
    return parser
def main_iSNP_per_patient(args):
    '''Aggregate tabular iSNP data per patient x position (all time points averaged)'''
    header = ['pos','patient','alleles','iSNV_freq','Hw',
        'eff_type','eff_codon_dna','eff_aa','eff_aa_pos','eff_prot_len','eff_gene','eff_protein']
    with open(args.outFile, 'wt') as outf:
        outf.write('\t'.join(header)+'\n')
        for row in iSNP_per_patient(util.file.read_tabfile_dict(args.inFile)):
            outf.write('\t'.join(map(str, [row.get(h,'') for h in header]))+'\n')
    return 0
__commands__.append(('iSNP_per_patient', parser_iSNP_per_patient))


#  ===================================================

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
