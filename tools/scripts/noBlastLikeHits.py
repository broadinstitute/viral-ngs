#!/usr/bin/env python
'''Written by Kristian Andersen.'''

import argparse
import sys

parser = argparse.ArgumentParser(description='This program outputs to stdout reads that have no blast hits')
parser.add_argument('-b', action="store", dest="blastPath", required=True, help="path to the blast-like hits file")
parser.add_argument('-r', action="store", dest="readsPath", required=True, help="path to the reads file")
parser.add_argument('-m',
                    action="store",
                    dest="hit",
                    required=True,
                    help="hit => output reads with hits, nohit => output reads with no hits",
                    choices=['hit', 'nohit'])
args = parser.parse_args()

# finds the nth and n+1th occurrence of a substring


def find_nth(str, substr, n):
    pos = list()
    i = 0
    for j in range(n + 1):
        i = str.find(substr, i + len(substr))
        if (j == n - 1):
            pos.append(i)
        elif (j == n):
            pos.append(i)
            break
    return pos


blastReads = {}
blastFile = open(args.blastPath, 'r')
for line in blastFile:
    if line.find('#') != 0:
        pos = find_nth(line, '\t', 6)
        blastReads[line[pos.pop(0) + 1:pos.pop(0)]] = True
blastFile.close()

FILE = open(args.readsPath, 'r')
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
else:
    print("Your file does not appear to be fasta or fastq")
    sys.exit()

if args.hit == "nohit":
    for i in range(0, len(offsets)):
        FILE.seek(offsets[i])
        line = FILE.readline()
        if not (line.split()[0][1:] in blastReads):
            if i == count - 1:
                FILE.seek(offsets[i])
                while True:
                    line = FILE.readline()
                    if not line:
                        break
                    sys.stdout.write(line)
            else:
                curOffset = 0
                targetOffset = offsets[i + 1] - offsets[i]
                FILE.seek(offsets[i])
                while curOffset != targetOffset:
                    line = FILE.readline()
                    curOffset += len(line)
                    sys.stdout.write(line)
else:
    for i in range(0, len(offsets)):
        FILE.seek(offsets[i])
        line = FILE.readline()
        if line.split()[0][1:] in blastReads:
            if i == count - 1:
                FILE.seek(offsets[i])
                while True:
                    line = FILE.readline()
                    if not line:
                        break
                    sys.stdout.write(line)
            else:
                curOffset = 0
                targetOffset = offsets[i + 1] - offsets[i]
                FILE.seek(offsets[i])
                while curOffset != targetOffset:
                    line = FILE.readline()
                    curOffset += len(line)
                    sys.stdout.write(line)
FILE.close()
