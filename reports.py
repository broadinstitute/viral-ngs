#!/usr/bin/env python
''' Reports
'''

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse, logging, subprocess, glob, os, os.path
import pysam
import util.cmd, util.file

log = logging.getLogger(__name__)



def consolidate_bamstats(inFiles, outFile):
    with util.file.open_or_gzopen(outFile, 'wt') as outf:
        header = []
        out_n = 0
        for fn in inFiles:
            out = {}
            with util.file.open_or_gzopen(fn, 'rt') as inf:
                for line in inf:
                    k,v = line.rstrip('\n').split('\t')
                    out[k] = v
                    if out_n==0:
                        header.append(k)
            if out_n==0:
                outf.write('\t'.join(header)+'\n')
            outf.write('\t'.join([out.get(h,'') for h in header])+'\n')
            out_n += 1

def parser_consolidate_bamstats() :
    parser = argparse.ArgumentParser(
        description='''Consolidate multiple bamstats reports into one.''')
    parser.add_argument('inFiles', help='Input report files.', nargs='+')
    parser.add_argument('outFile', help='Output report file.')
    return parser
def main_consolidate_bamstats(args) :
    consolidate_bamstats(args.inFiles, args.outFile)
    return 0
__commands__.append(('consolidate_bamstats',
    main_consolidate_bamstats, parser_consolidate_bamstats))



def consolidate_fastqc(inDirs, outFile):
    with util.file.open_or_gzopen(outFile, 'wt') as outf:
        header = ['Sample']
        out_n = 0
        for sdir in inDirs:
            out = {}
            with open(os.path.join(sdir, 'summary.txt'), 'rt') as inf:
                for line in inf:
                    v,k,fn = line.strip().split('\t')
                    out[k] = v
                    if out_n==0:
                        header.append(k)
                    if not fn.endswith('.bam'):
                        raise
                    out['Sample'] = fn[:-len('.bam')]
            if out_n==0:
                outf.write('\t'.join(header)+'\n')
            outf.write('\t'.join([out.get(h,'') for h in header])+'\n')
            out_n += 1

def parser_consolidate_fastqc() :
    parser = argparse.ArgumentParser(
        description='''Consolidate multiple FASTQC reports into one.''')
    parser.add_argument('inDirs', help='Input FASTQC directories.', nargs='+')
    parser.add_argument('outFile', help='Output report file.')
    return parser
def main_consolidate_fastqc(args) :
    consolidate_fastqc(args.inDirs, args.outFile)
    return 0
__commands__.append(('consolidate_fastqc',
    main_consolidate_fastqc, parser_consolidate_fastqc))



def consolidate_coverage(inFiles, adj, outFile):
    ending = '.coverage_%s.txt' % adj
    with util.file.open_or_gzopen(outFile, 'wt') as outf:
        for fn in inFiles:
            if not fn.endswith(ending):
                raise Exception()
            s = fn.split('/')[-1][:-len(ending)]
            with open(fn, 'rt') as inf:
                for line in inf:
                    outf.write(line.rstrip('\n') + '\t' + s + '\n')

def parser_consolidate_coverage() :
    parser = argparse.ArgumentParser(
        description='''Consolidate multiple coverage reports into one.''')
    parser.add_argument('inFiles', help='Input coverage files.', nargs='+')
    parser.add_argument('adj', help='Report adjective.')
    parser.add_argument('outFile', help='Output report file.')
    return parser
def main_consolidate_coverage(args) :
    consolidate_coverage(args.inFiles, args.adj, args.outFile)
    return 0
__commands__.append(('consolidate_coverage',
    main_consolidate_coverage, parser_consolidate_coverage))



def consolidate_spike_count(inFiles, outFile):
    with open(outFile, 'wt') as outf:
        for fn in inFiles:
            s = fn.split('/')[-1]
            if not s.endswith('.spike_count.txt'):
                raise Exception()
            s = s[:-len('.spike_count.txt')]
            with open(fn, 'rt') as inf:
                for line in inf:
                    if not line.startswith('Input bam'):
                        spike, count = line.strip().split('\t')
                        outf.write('\t'.join([s, spike, count])+'\n')

def parser_consolidate_spike_count() :
    parser = argparse.ArgumentParser(
        description='''Consolidate multiple spike count reports into one.''')
    parser.add_argument('inFiles', help='Input coverage files.', nargs='+')
    parser.add_argument('outFile', help='Output report file.')
    return parser
def main_consolidate_spike_count(args) :
    consolidate_spike_count(args.inFiles, args.outFile)
    return 0
__commands__.append(('consolidate_spike_count',
    main_consolidate_spike_count, parser_consolidate_spike_count))



# =======================

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
