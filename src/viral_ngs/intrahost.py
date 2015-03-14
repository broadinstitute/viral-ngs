#!/usr/bin/env python
'''This script contains a number of utilities for intrahost variant calling
and annotation for viral genomes.
'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org, "\
             "swohl@broadinstitute.org, irwin@broadinstitute.org"
__commands__ = []

import argparse, logging, itertools, re, shutil, tempfile, os
import Bio.AlignIO, Bio.SeqIO, Bio.Data.IUPACData
import pysam
import util.cmd, util.file, util.vcf, util.misc
from util.stats import mean, median, fisher_exact, chi2_contingency
from interhost import CoordMapper
from tools.vphaser2 import Vphaser2Tool
from tools.samtools import SamtoolsTool

log = logging.getLogger(__name__)

#  ============= class AlleleFieldParser =================

class AlleleFieldParser(object) :
    """
    Class for converting between the string and list representation
    of the fields in the allele columns of vphaser_one_sample output
    (corresponding to the SNP_or_LP_Profile columns in the V-Phaser 2 output).
    """
    def __init__(self, field = None,
                 allele = None, fcount = None, rcount = None, libBiasPval = None,
                 libCounts = None) :
        """ Input is either the string stored in one of the allele columns.
            or set field values.
        """
        if field == None :
            self._allele = allele
            self._strandCounts = [fcount, rcount]
            self._libBiasPval = libBiasPval
            self._libCounts = libCounts # libCounts is a list of 2-element lists
        else :
            words = field.split(':')
            self._allele = words[0]
            self._strandCounts = [int(words[1]), int(words[2])]
            self._libCounts = [[int(words[ii]), int(words[ii + 1])]
                              for ii in range(3, len(words) - 1, 2)]
            self._libBiasPval = float(words[-1])
   
    def __repr__(self) :
        """Convert to string representation."""
        return ':'.join([self._allele] +
                        list(map(str, self._strandCounts)) +
                        sum((list(map(str, libCount))
                             for libCount in self._libCounts),
                            []) +
                        ['%.4g' % self._libBiasPval])

    def allele(self) :
        """ Return allele:
                A, C, G, or T for SNVs,
                Dn with n > 0 for deletions,
                Ibases where bases is a string of two or more bases for inserts
        """
        return self._allele
    
    def total(self) :
        return sum(self._strandCounts)

    def strand_counts(self) :
        "Return [# forward reads (all libs), # reverse reads (all libs)]"
        return self._strandCounts
    
    def allele_and_strand_counts(self) :
        return [self.allele()] + self.strand_counts()

    def lib_counts(self) :
        """Yield [# forward reads, # reverse reads] for each library, in order
            of read groups in BAM file.
        """
        for counts in self._libCounts :
            yield counts

    def lib_bias_pval(self) :
        "Return a p-value on whether there is a library bias for this allele."
        return self._libBiasPval

#  ========== vphaser_one_sample =================

defaultMinReads = 5
defaultMaxBias = 10

def vphaser_one_sample(inBam, inConsFasta, outTab, vphaserNumThreads = None,
                       minReadsEach = None, maxBias = None) :
    ''' Input: a single BAM file, representing reads from one sample, mapped to
            its own consensus assembly. It may contain multiple read groups and 
            libraries.
        Output: a tab-separated file with no header containing filtered
            V Phaser-2 output variants with additional column for
            sequence/chrom name, and library counts and p-values appended to
            the counts for each allele.
    '''
    if minReadsEach != None and minReadsEach < 0:
        raise Exception('minReadsEach must be at least 0.')
    variantIter = Vphaser2Tool().iterate(inBam, vphaserNumThreads)
    filteredIter = filter_strand_bias(variantIter, minReadsEach, maxBias)
    libraryFilteredIter = compute_library_bias(filteredIter, inBam, inConsFasta)
    with util.file.open_or_gzopen(outTab, 'wt') as outf :
        for row in libraryFilteredIter :
            outf.write('\t'.join(row) + '\n')

def filter_strand_bias(isnvs, minReadsEach = None, maxBias = None) :
    ''' Take an iterator of V-Phaser output (plus chromosome name prepended)
        and perform hard filtering for strand bias
    '''
    alleleCol = 7 # First column of output with allele counts
    if minReadsEach == None :
        minReadsEach = defaultMinReads
    if maxBias == None :
        maxBias = defaultMaxBias
    for row in isnvs:
        front = row[:alleleCol]
        for fieldInd in range(len(row) - 1, alleleCol - 1, -1) :
            f, r = AlleleFieldParser(row[fieldInd]).strand_counts()
            if  (int(f)<minReadsEach or int(r)<minReadsEach
                    or (minReadsEach > 0 and not
                        (maxBias >= (float(f)/float(r)) >= 1.0/maxBias))) :
                del row[fieldInd]
        if len(row) > alleleCol + 1:
            row[alleleCol:] = sorted(row[alleleCol:],
                key = lambda field : AlleleFieldParser(field).total(),
                reverse = True)
            mac = sum(AlleleFieldParser(field).total() for field in row[alleleCol + 1:])
            tot = sum(AlleleFieldParser(field).total() for field in row[alleleCol:])
            row[2] = AlleleFieldParser(row[alleleCol + 1]).allele()
            row[3] = AlleleFieldParser(row[alleleCol]).allele()
            row[6] = '%.6g' % (100.0*mac/tot)
            yield row

def compute_library_bias(isnvs, inBam, inConsFasta) :
    ''' For each variant, compute read counts in each library and p-value for
          library bias; append them to string for each variant.
        Format is allele:totalF:totalR:1stLibFCount:1stLibRCount:2ndLibFCount:...:p-val.
        Library counts are in same order as in bam header.
        Note: Total was computed by vphaser, library counts by samtools mpileup,
          so total might not be sum of library counts.
    '''
    alleleCol = 7 # First column of output with allele counts
    samtoolsTool = SamtoolsTool()
    readGroups = list(r for r in samtoolsTool.getReadGroups(inBam)
        if samtoolsTool.count(inBam, ['-r', r])>0)
    tempDir = tempfile.mkdtemp()
    libBams = [os.path.join(tempDir, groupID + '.bam') for groupID in readGroups]
    for groupID, libBam in zip(readGroups, libBams) :
        samtoolsTool.view(['-r', groupID, '-b'], inBam, libBam)
        samtoolsTool.index(libBam)
    for row in isnvs :
        consensusAllele = row[3]
        pos = int(row[1]) if consensusAllele != 'i' else int(row[1]) - 1
        chrom = row[0]
        libCounts = [get_mpileup_allele_counts(libBam, chrom, pos, inConsFasta)
                     for libBam in libBams]
        numAlleles = len(row) - alleleCol
        countsMatrix = [[0] * numAlleles for lib in libBams]
        libCountsByAllele = []
        for alleleInd in range(numAlleles) :
            allele = row[alleleCol + alleleInd].split(':')[0]
            libCountsByAllele.append([])
            for libAlleleCounts, countsRow in zip(libCounts, countsMatrix) :
                f, r = libAlleleCounts.get(allele, [0, 0])
                libCountsByAllele[-1].append([f, r])
                countsRow[alleleInd] += f + r
        for alleleInd in range(numAlleles) :
            contingencyTable = [
                [         countsRow[alleleInd]         for countsRow in countsMatrix],
                [sum(countsRow) - countsRow[alleleInd] for countsRow in countsMatrix]]
            rowSums = map(sum, contingencyTable)
            dofs = len(libCounts) - 1
            if dofs < 1 :
                pval = 1.0
            elif min(rowSums) ** dofs / dofs < 10000 :
                # At this cutoff, fisher_exact should take <~ 0.1 sec
                pval = fisher_exact(contingencyTable)
            else :
                pval = chi2_contingency(contingencyTable)
            row[alleleCol + alleleInd] = str(AlleleFieldParser(None,
                *(row[alleleCol + alleleInd].split(':') +
                  [pval, libCountsByAllele[alleleInd]])))
        yield row
    shutil.rmtree(tempDir)

def parse_alleles_string(allelesStr) :
    # Return {allele : [forwardCount, reverseCount]}
    # For reference, allele is '.' rather than real allele
    alleleCounts = {} # allele : [forwardCount, reverseCount]
    pos = -1
    digits = re.compile('[0-9]+')
    while pos < len(allelesStr) - 1 :
        pos += 1
        c = allelesStr[pos]
        if c in '.,' :
            allele = '.'
            isRev = c == ','
        elif c in '<>$*' : # Reference skip, end of read, placeholder
            continue       # Not interested in these
        elif c == '^' : # Start of read
            pos += 1    # Skip quality character
            continue
        elif c in 'ACGTNacgtn' :
            allele = c.upper()
            isRev = c == c.lower()
        elif c in '+-' : # e.g., +3aaa
            mat = digits.match(allelesStr, pos + 1)
            indelLen = int(allelesStr[mat.start() : mat.end()])
            indelStr = allelesStr[mat.end() : mat.end() + indelLen]
            allele = 'I' + indelStr.upper() if c == '+' else 'D' + str(indelLen)
            isRev = indelStr == indelStr.lower()
            pos += mat.end() - mat.start() + indelLen
        else :
            raise Exception('Unknown allele type %s' % c)
        alleleCounts.setdefault(allele, [0, 0])
        alleleCounts[allele][isRev] += 1
    return alleleCounts

def get_mpileup_allele_counts(inBam, chrom, pos, inConsFasta) :
    """ Return {allele : [forwardCount, reverseCount], ...}
        allele is:
            Iins for insertions where ins represents the inserted bases
            Dlen for deletions where len is the length of the deletion
            base itself for non-indels
            'i' or 'd', in which case report count for consensus.
    """
    pileupFileName = util.file.mkstempfname('.txt')
    SamtoolsTool().mpileup(inBam, pileupFileName,
                           ['-A',
                            '-r', '%s:%d-%d' % (chrom, pos, pos),
                            '-B', '-d', '50000', '-L', '50000', '-Q', '0',
                            '-f', inConsFasta])
    with open(pileupFileName) as pileupFile :
        words = pileupFile.readline().split('\t')
    if len(words)<5:
        # empty output files means no reads pile up on this position
        return {}
    alleleCounts = parse_alleles_string(words[4])

    # '.' is whatever mpileup thinks is the reference base (which might be
    # different from vphaser's consensus base). This is the count we want to
    # report for this base, but we also want to report this count for what
    # vphaser calls 'i' or 'd'.
    # This will probably be wrong in the unlikely case that there are both
    # an indel and a snp at the same position, and vphaser's consensus base
    # is different from mpileup's at that position.
    refAllele = words[2]
    alleleCounts['i'] = alleleCounts['d'] = alleleCounts[refAllele] = \
        alleleCounts.get('.', [0, 0])
    
    return alleleCounts

def parser_vphaser_one_sample(parser = argparse.ArgumentParser()) :
    parser.add_argument("inBam",
        help = "Input Bam file.")
    parser.add_argument("inConsFasta",
        help = "Consensus assembly fasta.")
    parser.add_argument("outTab", help = "Tab-separated headerless output file.")
    parser.add_argument("--vphaserNumThreads", type = int, default = None,
        help="Number of threads in call to V-Phaser 2.")
    parser.add_argument("--minReadsEach", type = int, default = defaultMinReads,
        help = "Minimum number of reads on each strand (default: %(default)s).")
    parser.add_argument("--maxBias", type = int, default = defaultMaxBias,
        help = """Maximum allowable ratio of number of reads on the two strands
                (default: %(default)s). Ignored if minReadsEach = 0.""")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, vphaser_one_sample, split_args = True)
    return parser
__commands__.append(('vphaser_one_sample', parser_vphaser_one_sample))

#  ========== vphaser =================

def parser_vphaser(parser = argparse.ArgumentParser()) :
    parser.add_argument("inBam",
        help = "Input Bam file.")
    parser.add_argument("outTab", help = "Tab-separated headerless output file.")
    parser.add_argument("--numThreads", type = int, default = None,
        help="Number of threads in call to V-Phaser 2.")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, vphaser_main, split_args = True)
    return parser
def vphaser_main(inBam, outTab, numThreads = None) :
    """ Run V-Phaser 2 on the input file without any additional filtering.
        Combine the non-header lines of the CHROM.var.raw.txt files it produces, 
            adding CHROM as the first field on each line.
    """
    with open(outTab, 'wt') as outf :
        for row in Vphaser2Tool().iterate(inBam, numThreads) :
            outf.write('\t'.join(row) + '\n')
__commands__.append(('vphaser', parser_vphaser))

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

def strip_accession_version(acc):
    ''' If this is a Genbank accession with a version number,
        remove the version number.
    '''
    m = re.match(r"^(\S+)\.\d+$", acc)
    if m:
        acc = m.group(1)
    return acc

def merge_to_vcf(refFasta, outVcf, samples, isnvs, assemblies, strip_chr_version=False, naive_filter=False):
    ''' Combine and convert vPhaser2 parsed filtered output text files into VCF format.
        Assumption: consensus assemblies do not extend beyond ends of reference.
    '''
    
    # setup
    if not (len(samples) == len(isnvs) == len(assemblies)):
        raise LookupError("samples, isnvs, and assemblies must have the same number of elements")
    samp_to_seqIndex = dict((s, Bio.SeqIO.index(fasta, 'fasta'))
        for s, fasta in zip(samples, assemblies))
    samp_to_isnv = dict(zip(samples, isnvs))
    samp_to_cmap = dict((s, CoordMapper(genome, refFasta))
        for s, genome in zip(samples, assemblies))
    with open(refFasta, 'rU') as inf:
        ref_chrlens = list((seq.id, len(seq)) for seq in Bio.SeqIO.parse(inf, 'fasta'))
    
    if outVcf.endswith('.vcf.gz'):
        tmpVcf = util.file.mkstempfname('.vcf')
    elif outVcf.endswith('.vcf'):
        tmpVcf = outVcf
    else:
        raise ValueError("outVcf must end in .vcf or .vcf.gz")
    
    log.info("loaded CoordMapper for all genomes, starting VCF merge...")
    with open(tmpVcf, 'wt') as outf:
    
        # write header
        outf.write('##fileformat=VCFv4.1\n')
        outf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outf.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        outf.write('##FORMAT=<ID=NL,Number=R,Type=Integer,Description="Number of libraries observed per allele">\n')
        outf.write('##FORMAT=<ID=LB,Number=R,Type=Float,Description="Library bias observed per allele (Fishers Exact P-value)">\n')
        for c, clen in ref_chrlens:
            if strip_chr_version:
                c = strip_accession_version(c)
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
                        s_chrom = samp_to_cmap[s].mapBtoA(ref_sequence.id)
                        if row[0] == s_chrom:
                            allele_fields = list(AlleleFieldParser(x) for x in row[7:] if x)
                            row = {'sample':s, 'CHROM':ref_sequence.id,
                                's_chrom':s_chrom, 's_pos':int(row[1]),
                                's_alt':row[2], 's_ref': row[3],
                                'alleles':list(x.allele_and_strand_counts() for x in allele_fields),
                                'n_libs':dict(
                                    (x.allele(), sum(1 for f,r in x.lib_counts() if f+r>0))
                                    for x in allele_fields),
                                'lib_bias':dict(
                                    (x.allele(), x.lib_bias_pval())
                                    for x in allele_fields),
                            }
                            # make a sorted allele list
                            row['allele_counts'] = list(sorted(
                                [(a,int(f)+int(r)) for a,f,r in row['alleles']],
                                key=(lambda x:x[1]), reverse=True))
                            # naive filter (quick and dirty)
                            if naive_filter:
                                tot_n = sum(n for a,n in row['allele_counts'])
                                row['allele_counts'] = list((a,n)
                                    for a,n in row['allele_counts']
                                    if row['n_libs'][a]>=2 and float(n)/tot_n >= 0.005)
                            # reposition vphaser deletions minus one to be consistent with
                            # VCF conventions
                            if row['s_alt'].startswith('D'):
                                for a,n in row['allele_counts']:
                                    if a[0] not in ('D','i'):
                                        log.error("allele_counts: " + str(row['allele_counts']))
                                        raise Exception("deletion alleles must always start with D or i")
                                row['s_pos'] = row['s_pos']-1
                            # map position back to reference coordinates
                            row['POS'] = samp_to_cmap[s].mapAtoB(s_chrom, row['s_pos'], side = -1)[1]
                            row['END'] = samp_to_cmap[s].mapAtoB(s_chrom, row['s_pos'], side = 1)[1]
                            if row['POS'] == None or row['END'] == None:
                                raise Exception('consensus extends beyond start or end of reference.')
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
                        end = max(end, row['END'])
                        for a,n in row['allele_counts']:
                            if a.startswith('D'):
                                # end of deletion in sample's coord space
                                local_end = row['s_pos']+int(a[1:])
                                # end of deletion in reference coord space
                                ref_end = samp_to_cmap[row['sample']].mapAtoB(row['s_chrom'], local_end, side = 1)[1]
                                if ref_end == None:
                                    raise Exception('consensus extends ' \
                                     'beyond start or end of reference.')
                                end = max(end, ref_end)
                
                    # find reference allele and consensus alleles
                    refAllele = str(ref_sequence[pos-1:end].seq)
                    consAlleles = {} # the full pos-to-end consensus assembly sequence for each sample
                    samp_offsets = {} # {sample : isnv's index in its consAllele string}
                    for row in rows :
                        s_pos = row['s_pos']
                        sample = row['sample']
                        if samp_offsets.get(sample, s_pos) != s_pos :
                            raise NotImplementedError('Sample %s has variants at 2 '
                                'positions %s mapped to same reference position (%s)' %
                                (sample, (s_pos, samp_offsets[sample]), pos))
                        samp_offsets[sample] = s_pos
                    for s in samples:
                        cons_start = samp_to_cmap[s].mapBtoA(ref_sequence.id, pos, side = -1)[1]
                        cons_stop  = samp_to_cmap[s].mapBtoA(ref_sequence.id, end, side =  1)[1]
                        if cons_start == None or cons_stop == None :
                            log.warning("dropping consensus because allele is outside "
                                "consensus for %s at %s-%s." % (s, pos, end))
                            continue
                        cons = samp_to_seqIndex[s][samp_to_cmap[s].mapBtoA(ref_sequence.id)]
                        allele = str(cons[cons_start-1:cons_stop].seq).upper()
                        if s in samp_offsets:
                            samp_offsets[s] -= cons_start
                        if all(a in set(('A','C','T','G')) for a in allele):
                            consAlleles[s] = allele
                        else:
                            log.warning("dropping unclean consensus for %s at %s-%s: %s" % (s, pos, end, allele))
                    
                    # define genotypes and fractions
                    iSNVs = {} # {sample : {allele : fraction, ...}, ...}
                    iSNVs_n_libs = {} # {sample : {allele : n libraries > 0, ...}, ...}
                    iSNVs_lib_bias = {} # {sample : {allele : pval, ...}, ...}
                    for s in samples:
                        
                        # get all rows for this sample and merge allele counts together
                        acounts = dict(itertools.chain.from_iterable(row['allele_counts']
                            for row in rows if row['sample'] == s))
                        nlibs = dict(itertools.chain.from_iterable(row['n_libs'].items()
                            for row in rows if row['sample'] == s))
                        libbias = dict(itertools.chain.from_iterable(row['lib_bias'].items()
                            for row in rows if row['sample'] == s))
                        if 'i' in acounts and 'd' in acounts:
                            # This sample has both an insertion line and a deletion line at the same spot!
                            # To keep the reference allele from be counted twice, once as an i and once
                            # as a d, average the counts and get rid of one of them.
                            acounts['i'] = int(round((acounts['i'] + acounts['d'])/2.0,0))
                            del acounts['d']
                            nlibs['i'] = max(nlibs['i'], nlibs['d'])
                            libbias['i'] = max(libbias['i'], libbias['d'])
                        
                        if acounts and s in consAlleles:
                            # we have iSNV data on this sample
                            consAllele = consAlleles[s]
                            tot_n = sum(acounts.values())
                            iSNVs[s] = {} # {allele : fraction, ...}
                            iSNVs_n_libs[s] = {}
                            iSNVs_lib_bias[s] = {}
                            for orig_a,n in acounts.items():
                                f = float(n)/tot_n
                                a = orig_a
                                if a.startswith('I'):
                                    # insertion point is relative to each sample
                                    insert_point = samp_offsets[s]+1
                                    a = consAllele[:insert_point] + a[1:] + consAllele[insert_point:]
                                elif a.startswith('D'):
                                    # deletion is the first consensus base, plus remaining
                                    # consensus seq with the first few positions dropped off
                                    cut_left = samp_offsets[s]+1
                                    cut_right = samp_offsets[s]+1+int(a[1:])
                                    a = consAllele[:cut_left] + consAllele[cut_right:]
                                elif a in ('i','d'):
                                    # this is vphaser's way of saying the "reference" (majority/consensus)
                                    # allele, in the face of other indel variants
                                    a = consAllele
                                else:
                                    # this is a SNP
                                    if a not in set(('A','C','T','G')):
                                        raise Exception()
                                    if f>0.5 and a!=consAllele[samp_offsets[s]]:
                                        log.warning("vPhaser and assembly pipelines mismatch at "
                                            "%s:%d %s - consensus %s, vPhaser %s" %
                                            (ref_sequence.id, pos, s, consAllele[0], a))
                                    new_allele = list(consAllele)
                                    new_allele[samp_offsets[s]] = a
                                    a = ''.join(new_allele)
                                if not (a and a==a.upper()):
                                    raise Exception()
                                iSNVs[s][a] = f
                                iSNVs_n_libs[s][a] = nlibs[orig_a]
                                iSNVs_lib_bias[s][a] = libbias[orig_a]
                            if all(len(a)==1 for a in iSNVs[s].keys()):
                                if consAllele not in iSNVs[s]:
                                    raise Exception()
                        elif s in consAlleles:
                            # there is no iSNV data for this sample, so substitute the consensus allele
                            iSNVs[s] = {consAlleles[s]:1.0}

                    # get unique alleles list for this position, in this order:
                    # first:   reference allele,
                    # next:    consensus allele for each sample, in descending order of
                    #          number of samples with that consensus,
                    # finally: all other alleles, sorted first by number of containing samples,
                    #          then by intrahost read frequency summed over the population,
                    #          then by the allele string itself.
                    alleles_cons = [a
                                    for a,n in sorted(util.misc.histogram(consAlleles.values()).items(),
                                                      key=lambda x:x[1], reverse=True)
                                    if a!=refAllele]
                    alleles_isnv = list(itertools.chain.from_iterable([iSNVs[s].items() for s in samples if s in iSNVs]))
                    alleles_isnv2 = []
                    for a in set(a for a,n in alleles_isnv):
                        counts = list(x[1] for x in alleles_isnv if x[0]==a)
                        if len(counts)>0 and sum(counts)>0:
                            # if we filtered any alleles above, make sure to omit absent alleles
                            alleles_isnv2.append((len(counts),sum(counts),a))
                        else:
                            log.info("dropped allele {} at position {}:{}".format(a, ref_sequence.id, pos))
                    alleles_isnv = list(allele for n_samples, n_reads, allele in reversed(sorted(alleles_isnv2)))
                    alleles = list(util.misc.unique([refAllele] + alleles_cons + alleles_isnv))
                    
                    # map alleles from strings to numeric indexes
                    if not alleles:
                        raise Exception()
                    elif len(alleles)==1:
                        # if we filtered any alleles above, skip this position if there is no variation left here
                        log.info("dropped position {}:{} due to lack of variation".format(ref_sequence.id, pos))
                        continue
                    alleleMap = dict((a,i) for i,a in enumerate(alleles))
                    genos = [str(alleleMap.get(consAlleles.get(s),'.')) for s in samples]
                    freqs = [(s in iSNVs) and ','.join(map(str, [iSNVs[s].get(a,0.0) for a in alleles[1:]])) or '.'
                             for s in samples]
                    nlibs = [(s in iSNVs_n_libs) and ','.join([str(iSNVs_n_libs[s].get(a,0)) for a in alleles]) or '.'
                             for s in samples]
                    pvals = [(s in iSNVs_lib_bias) and ','.join([str(iSNVs_lib_bias[s].get(a,'.')) for a in alleles]) or '.'
                             for s in samples]

                    # prepare output row and write to file
                    c = ref_sequence.id
                    if strip_chr_version:
                        c = strip_accession_version(c)
                    out = [c, pos, '.',
                        alleles[0], ','.join(alleles[1:]),
                        '.', '.', '.', 'GT:AF:NL:LB']
                    out = out + list(map(':'.join, zip(genos, freqs, nlibs, pvals)))
                    outf.write('\t'.join(map(str, out))+'\n')
    
    # close fasta files
    for idx in samp_to_seqIndex.values():
        idx.close()
    
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
    parser.add_argument("--strip_chr_version",
        default=False, action="store_true", dest="strip_chr_version",
        help="""If set, strip any trailing version numbers from the
            chromosome names. If the chromosome name ends with a period
            followed by integers, this is interepreted as a version number
            to be removed. This is because Genbank accession numbers are
            often used by SnpEff databases downstream, but without the
            corresponding version number.  Default is false (leave
            chromosome names untouched).""")
    parser.add_argument("--naive_filter",
        default=False, action="store_true", dest="naive_filter",
        help="""If set, keep only the alleles that have at least
            two independent libraries of support and allele freq > 0.005.
            Default is false (do not filter at this stage).""")
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


#  =================== iSNV_table ================================

def parse_eff(eff_field):
    ''' parse the old snpEff "EFF" INFO field '''
    out = {}
    effs = [eff.rstrip(')').replace('(','|').split('|') for eff in eff_field.split(',')]
    effs = [[eff[i] for i in (0,3,4,5,6,9,11)] for eff in effs]
    effs = [eff for eff in effs if eff[5] not in ('sGP','ssGP') and int(eff[6])<2]
    if not len(effs)==1:
        raise Exception()
    eff = effs[0]
    if eff[2]:
        aa = eff[2].split('/')[0]
        assert aa.startswith('p.')
        aa = aa[2:]
        m = re.search(r"(\d+)", aa)
        out['eff_aa_pos'] = int(m.group(1))
    (out['eff_type'], out['eff_codon_dna'], out['eff_aa'], out['eff_prot_len'], out['eff_gene'], out['eff_protein'], rank) = eff
    return out

def parse_ann(ann_field, alleles, transcript_blacklist=set(('GP.2','GP.3'))):
    ''' parse the new snpEff "ANN" INFO field '''
    
    if not all(len(a)==1 for a in alleles):
        ''' This is incredibly annoying: instead of talking about indel alleles using
            the VCF allele definitions or in the VCF order, the ANN field recomputes the
            part of the allele that looks distinct from the reference and talks about
            that instead. Note that this means some alleles (particularly deletions)
            will collapse and become non-unique.
            Oh, and py2 and py3 have different names for the same function, just for fun.
        '''
        ziplong = itertools.zip_longest if 'zip_longest' in dir(itertools) else itertools.izip_longest
        outalleles = []
        for a in alleles[1:]:
            if len(a)==len(alleles[0]):
                out_a = []
                for ref, alt in reversed(list(zip(alleles[0], a))):
                    if ref!=alt or out_a:
                        out_a.append(alt)
                a = ''.join(reversed(out_a))
            out_a = []
            for ref, alt in ziplong(alleles[0], a, fillvalue=''):
                if ref!=alt or out_a:
                    out_a.append(alt)
            outalleles.append(''.join(out_a))
        alleles = outalleles
    else:
        alleles = alleles[1:]
    
    effs = [eff.split('|') for eff in ann_field.split(',')]
    effs = [(eff[0], dict((k,eff[i]) for k,i in
            (('eff_type',1),('eff_gene',3),('eff_protein',6),
            ('eff_codon_dna',9),('eff_aa',10),
            ('eff_aa_pos',13),('eff_prot_len',13))))
        for eff in effs
        if eff[6] not in transcript_blacklist]
    effs_dict = dict(effs)
    if not effs:
        return {}
    if len(effs) != len(effs_dict):
        raise Exception("ANN field has non-unique alleles")
    for a in alleles:
        if a not in effs_dict:
           raise Exception("ANN field is missing ALT allele: " + a)
    if len(effs) != len(set(alleles)):
        raise Exception("ANN field has %s entries, but ALT field has %s unique alleles: %s" % (
            len(effs), len(set(alleles)), ','.join(alleles)))
    
    out = {}
    for k in ('eff_type', 'eff_codon_dna', 'eff_aa', 'eff_aa_pos', 'eff_prot_len', 'eff_gene', 'eff_protein'):
        a_out = []
        for a in alleles:
            v = effs_dict[a][k]
            if k=='eff_codon_dna' and v.startswith('c.'):
                v = v[2:]
            elif k=='eff_aa' and v.startswith('p.'):
                v = v[2:]
            elif k=='eff_aa_pos' and '/' in v:
                v = v.split('/')[0]
            elif k=='eff_prot_len' and '/' in v:
                v = v.split('/')[1]
            elif k=='eff_protein' and v=='GP.1':
                v = 'Glycoprotein'
            if v:
                a_out.append(v)
        out[k] = ','.join(util.misc.unique(a_out))
    return out

def iSNV_table(vcf_iter):
    for row in vcf_iter:
        info = dict(kv.split('=') for kv in row['INFO'].split(';') if kv and kv != '.')
        samples = [k for k in row.keys() if k not in set(('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'))]
        try:
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
                        for k,v in parse_eff(info['EFF']).items():
                            out[k] = v
                    if 'ANN' in info:
                        for k,v in parse_ann(info['ANN'], alleles=out['alleles'].split(',')).items():
                            out[k] = v
                    if 'PI' in info:
                        out['Hs_snp'] = info['PI']
                    if 'FWS' in info:
                        out['Fws_snp'] = info['FWS']
                    yield out
        except:
            log.error("VCF parsing error at {}:{}".format(row['CHROM'], row['POS']))
            raise

def parser_iSNV_table(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outFile", help="Output text file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_iSNV_table)
    return parser
def main_iSNV_table(args):
    '''Convert VCF iSNV data to tabular text'''
    header = ['chr','pos','sample','patient','time','alleles','iSNV_freq','Hw',
        'eff_type','eff_codon_dna','eff_aa','eff_aa_pos','eff_prot_len','eff_gene','eff_protein']
    with util.file.open_or_gzopen(args.outFile, 'wt') as outf:
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
