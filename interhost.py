#!/usr/bin/env python
''' This script contains a number of utilities for SNP calling, multi-alignment,
    phylogenetics, etc.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import Bio.AlignIO
import argparse, logging
import tools.muscle
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)


class CoordMapper:
    ''' Maps coordinates between genome A and genome B.
        The two genomes are described by fasta files and
        may contain multiple chromosomes, but must have the
        same number of chromosomes and must be described in
        the same order.
    '''
    def __init__(self, fastaA, fastaB, aligner=tools.muscle.MuscleTool):
        self.fastaA = fastaA
        self.fastaB = fastaB
        self.aligner = aligner()
    def _align(self):
        raise NotImplementedError()
        # assert that fastaA and fastaB have the same non-zero number of sequences
        # save a list of chr names for each genome (the names may be different)
        # assert that they are unique within each genome
        # save a dict mapping chr names from A -> B and another for B -> A
        # for each chr number, extract that chr from both genomes into a standalone
        #    temp file fasta, align with self.aligner, then get the gapped/aligned
        #    fasta from the output of the aligner, and construct a gap-string
        #    (CIGAR or whatever) between the two genomes.
        # save a dict mapping chr names from A to gapstrings that map A -> B
        # and vice versa
        # delete temp files
    def mapAtoB(self, c, p):
        raise NotImplementedError()
        # map chr name from A -> B
        # retrieve A->B gapstring
        # map pos from A->B
        # return (new_c, new_p)
    def mapBtoA(self, c, p):
        raise NotImplementedError()
        


# modified version of rachel's call_snps_3.py follows
def call_snps_3(inFasta, outVcf, REF="KJ660346.2"):
    a=Bio.AlignIO.read(inFasta, "fasta")
    ref_idx = find_ref(a, REF)
    with open(outVcf, 'wt') as outf:
        outf.write(vcf_header(a))
        for row in make_vcf(a, ref_idx, REF):
            outf.write('\t'.join(map(str, row))+'\n')
def find_ref(a, ref):
    for i in range(len(a)):
        if a[i].id == ref:
            return i
    return -1
def vcf_header(a):
    header = "##fileformat=VCFv4.1\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##contig=<ID=\"KM034562\",length=18957>\n"
    header += '#' + '\t'.join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [x.id for x in a]) + '\n'
    return header
def make_vcf(a, ref_idx, chrom):
    bases=set(["A", "C", "G", "T"])
    for i in range(len(a[0])):
        alt = []
        for j in range(len(a)):
            if (a[j][i] != a[ref_idx][i]) and ((a[ref_idx][i] in bases) and (a[j][i] in bases)) and a[j][i] not in alt:
                alt.append(a[j][i])
        if len(alt) > 0:
            row = [chrom, i+1, '.', a[ref_idx][i], ','.join(alt), '.', '.', '.', 'GT']
            vars = []
            for k in range(len(a)):
                if a[k][i] == a[ref_idx][i]:
                    vars.append(0)
                elif a[k][i] not in bases:
                    vars.append(".")
                else:
                    for m in range(0, len(alt)):
                        if a[k][i] == alt[m]:
                            vars.append(m+1)
            yield row+vars
            

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
