#!/usr/bin/env python
''' This script contains a number of utilities for SNP calling, multi-alignment,
    phylogenetics, etc.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import Bio.AlignIO
import argparse, logging
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)





# rachel's call_snps_3.py follows
def call_snps_3(inFasta):
    REF = "EBOV_2014_G3686"
    CHROM = "KM034562"
    a=Bio.AlignIO.read(inFasta, "fasta")
    ref_idx = find_ref(a, REF)
    print_header(a)
    make_vcf(a, ref_idx, CHROM)
def print_header(a):
    print("##fileformat=VCFv4.1")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    print("##contig=<ID=\"KM034562\",length=18957>")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", end="")
    for i in range(0, len(a) - 1):
        print(a[i].id + "\t", end="")
    print(a[len(a) - 1].id + "\n", end="")
def find_ref(a, ref):
    for i in range(0, len(a)):
        if a[i].id == ref:
            return i
    return -1
def make_vcf(a, ref_idx, chrom):
    bases=set(["A", "C", "G", "T"])
    for i in range(0,len(a[0])):
        alt = []
        for j in range(0,len(a)):
            if (a[j][i] != a[ref_idx][i]) and ((a[ref_idx][i] in bases) and (a[j][i] in bases)) and a[j][i] not in alt:
                alt.append(a[j][i])
        if len(alt) > 0:
            print(chrom + "\t" + str(i+1) + "\t.\t" + a[ref_idx][i] + "\t", end="")
            print(",".join(alt), end = "")
            print("\t.\t.\t.\tGT\t", end = "")
            vars = []
            for k in range(0,len(a)):
                if a[k][i] == a[ref_idx][i]:
                    vars.append("0")
                elif a[k][i] not in bases:
                    vars.append(".")
                else:
                    for m in range(0, len(alt)):
                        if a[k][i] == alt[m]:
                            vars.append(str(m+1))
            print("\t".join(vars))


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
