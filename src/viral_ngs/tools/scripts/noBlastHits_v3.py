#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description=\
  'This program outputs to stdout reads that have no blast hits')
parser.add_argument('-b',action="store",dest="blastPath",required=True,\
  help="path to the blast hits file")
parser.add_argument('-r',action="store",dest="readsPath",required=True,\
  help="path to the reads file")
parser.add_argument('-m',action="store",dest="hit",required=True,\
  help="hit => output reads with hits, nohit => output reads with no hits",choices=['hit','nohit'])
args = parser.parse_args()

blastReads = {}
blastFile = open(args.blastPath, 'r')
for line in blastFile:
  blastReads[(line[0:line.find('\t')])] = True
blastFile.close()

readsFile = open(args.readsPath, 'r')
nohit =  args.hit == "nohit"
isFastq = args.readsPath.endswith( '.fastq' )
while True:
  line1 = readsFile.readline()
  line2 = readsFile.readline()
  if not line2:
    break
  line3 = ''
  line4 = ''
  if isFastq:
    line3 = readsFile.readline()
    if not line3:
      break
    line4 = readsFile.readline()
    if not line4:
      break
  if nohit != (line1[1:line1.find('\n')] in blastReads):
    sys.stdout.write(line1+line2+line3+line4)
readsFile.close()
