# !/usr/bin/env python
#
# filter_short_seqs.py - filter sequences below a minimum length
#
# requires python >= 2.7

import sys, argparse, Bio.SeqIO


def main():
    parser = parser_filter_short_seqs()
    args = parser.parse_args()
    return main_filter_short_seqs(args)

def parser_filter_short_seqs():
    parser = argparse.ArgumentParser(description = "Check a contig.fasta file, if it falls below <length> then output it to bad_out <location/name>; otherwise, output it to good_out <location/name>")
    parser.add_argument("input", help="input sequence file")
    parser.add_argument("length", help="minimum length for contig", type=int)
    parser.add_argument("bad_out", help="output file if contig is too short")
    parser.add_argument("good_out", help="output file if contig is ok")
    parser.add_argument("-f", "--format", help="Format for input sequence (default: fasta)", default="fasta")
    parser.add_argument("-of", "--output-format",
                        help="Format for output sequence (default: fasta)", default="fasta")
    return parser

def main_filter_short_seqs(args):
    s = Bio.SeqIO.read(args.input, args.format)
    if (len(s) < args.length):
        f = open(args.bad_out, "w")
        Bio.SeqIO.write(s, f, args.output_format)
        f.close()
    else:
        f = open(args.good_out, "w")
        Bio.SeqIO.write(s, f, args.output_format)
        f.close()

if __name__ == '__main__':
    main()
