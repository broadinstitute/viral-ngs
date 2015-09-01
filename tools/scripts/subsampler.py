#!/usr/bin/env python
import argparse
import random
import sys
import subprocess

parser = argparse.ArgumentParser(
    description='This program outputs to stdout a user-defined number of reads from given reads file[s]')
parser.add_argument('-mode', action="store", dest="mode", required=True, help="s => single-end, p => paired-end")
# parser.add_argument('-format',action="store",dest="format",required=True,\
#  help="fasta => fasta, fastq => fastq",choices=['fasta','fastq'])
parser.add_argument('-n',
                    action="store",
                    type=int,
                    dest="numSeq",
                    required=True,
                    help="specifies number of reads to output")
parser.add_argument('-in',
                    action="store",
                    dest="inputs",
                    nargs='+',
                    required=True,
                    help="specifies input files (put both input files separated by a space for paired end)")
parser.add_argument('-out',
                    action="store",
                    dest="outputs",
                    nargs='*',
                    required=False,
                    default=['/dev/stdout'],
                    help="specifies output files (default for single end is stdout)")
args = parser.parse_args()
#format = args.format
numSeq = args.numSeq
mode = args.mode
reads = args.inputs
outs = args.outputs
if mode == 's':
    if len(reads) != 1 or len(outs) > 1:
        print("please specify only one read file and at most one output file (if paired end, specify with -pe)")
        sys.exit()
elif mode == 'p':
    if len(reads) != 2 or len(outs) != 2:
        print("please specify both read files and two output files")
        sys.exit()
fileName1 = reads[0]
if fileName1.find(".gz") != -1:
    subprocess.call("gunzip " + fileName1, shell=True)
    fileName1 = fileName1[0:-3]
FILE = open(fileName1, 'r')
if mode == 'p':
    fileName2 = reads[1]
    if fileName2.find(".gz") != -1:
        subprocess.call("gunzip " + fileName2, shell=True)
        fileName2 = fileName1[0:-3]
    FILE2 = open(fileName2, 'r')
random.seed()
# name of the input file (fasta or fastq)
# assumes input file is standard fasta/fastq format
"""
calculate number of total reads
pseudorandomly determine which reads will be extracted
"""
count = 0
char = FILE.readline()[0:1]
FILE.seek(0)
offsets = list()
# fasta format
if char == '>':
    while True:
        line = FILE.readline()
        if not line:
            break
        if line[0:1] == char:
            count += 1
            offsets.append(FILE.tell() - len(line))
    if mode == 'p':
        offsets2 = list()
        while True:
            line = FILE2.readline()
            if not line:
                break
            if line[0:1] == char:
                offsets2.append(FILE2.tell() - len(line))
# fastq format
elif char == '@':
    while True:
        offsets.append(FILE.tell())
        line = FILE.readline()
        if not line:
            offsets.pop(count)
            break
        count += 1
        linesToPlus = 0
        while FILE.readline()[0:1] != '+':
            linesToPlus += 1
        for i in range(0, linesToPlus):
            FILE.readline()
    if mode == 'p':
        offsets2 = list()
        while True:
            offsets2.append(FILE2.tell())
            line = FILE2.readline()
            if not line:
                offsets2.pop(count)
                break
            linesToPlus = 0
            while FILE2.readline()[0:1] != '+':
                linesToPlus += 1
            for i in range(0, linesToPlus):
                FILE2.readline()
else:
    print("Your file does not appear to be a valid format")
    sys.exit()

if count < numSeq:
    numSeq = count
selected = list(range(0, count))
random.shuffle(selected)
selected = sorted(selected[0:numSeq])

if mode == 's':
    out = open(outs[0], 'w')
    for i in range(0, numSeq):
        if selected[i] == count - 1:
            FILE.seek(offsets[selected[i]])
            while True:
                line = FILE.readline()
                if not line:
                    break
                out.write(line)
        else:
            curOffset = 0
            targetOffset = offsets[selected[i] + 1] - offsets[selected[i]]
            FILE.seek(offsets[selected[i]])
            while curOffset != targetOffset:
                line = FILE.readline()
                curOffset += len(line)
                out.write(line)
    FILE.close()
    out.close()
    sys.exit()
# to avoid incessant checking
elif mode == 'p':
    out1 = open(outs[0], 'w')
    out2 = open(outs[1], 'w')
    for i in range(0, numSeq):
        if selected[i] == count - 1:
            FILE.seek(offsets[selected[i]])
            FILE2.seek(offsets2[selected[i]])
            while True:
                line = FILE.readline()
                line2 = FILE2.readline()
                if not line:
                    break
                out1.write(line)
                out2.write(line2)
        else:
            curOffset = 0
            targetOffset = offsets[selected[i] + 1] - offsets[selected[i]]
            FILE.seek(offsets[selected[i]])
            FILE2.seek(offsets2[selected[i]])
            while curOffset != targetOffset:
                line = FILE.readline()
                line2 = FILE2.readline()
                curOffset += len(line)
                out1.write(line)
                out2.write(line2)
    FILE.close()
    FILE2.close()
    out1.close()
    out2.close()
    sys.exit()
