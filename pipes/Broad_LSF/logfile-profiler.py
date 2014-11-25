#!/usr/bin/env python
''' This script contains utilities for analyzing LSF reports for
    data on runtime and performance.
'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, re, time, glob
import util.cmd, util.file

log = logging.getLogger(__name__)



def read_lsf_logfile(infname):
    out = {}
    num_dash_lines = 0
    with open(infname, 'rt') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('Job <'):
                mo = re.match(r'Job <(\w+?)-(\w+)> was submitted from host', line)
                out['job_rule'] = mo.group(1)
                out['job_suffix'] = mo.group(2)
            elif line.startswith('Job was executed on host'):
                mo = re.match(r'Job was executed on host(s) <(\w+)>, in queue <(\w+)>', line)
                out['exec_host'] = mo.group(1)
                out['queue'] = mo.group(2)
            elif line.startswith('Started at'):
                mo = re.match(r'Started at (.+)$', line)
                out['start_time'] = mo.group(1)
            elif line.startswith('Results reported at'):
                mo = re.match(r'Results reported at (.+)$', line)
                out['end_time'] = mo.group(1)
            elif line.startswith('----------------'):
                num_dash_lines += 1
            elif line.startswith('The output (if any) follows'):
                break
            else:
                if line and num_dash_lines == 2:
                    if 'status' not in out:
                        out['status'] = line
                    elif ':' in line:
                        k,v = [s.strip() for s in line.split(':')]
                        out[k]=v
    if 'start_time' in out and 'end_time' in out:
        duration = time.strptime(out['end_time']) - time.strptime(out['start_time'])
        out['run_time'] = duration.total_seconds()
    return out

def read_all_logfiles(dirname):
    for fname in glob.glob(dirname):
        yield read_lsf_logfile(fname)


def parser_filter_short_seqs():
    parser = argparse.ArgumentParser(description = "Check sequences in inFile, retaining only those that are at least minLength")
    parser.add_argument("inFile", help="input sequence file")
    parser.add_argument("minLength", help="minimum length for contig", type=int)
    parser.add_argument("outFile", help="output file")
    parser.add_argument("-f", "--format", help="Format for input sequence (default: fasta)", default="fasta")
    parser.add_argument("-of", "--output-format",
                        help="Format for output sequence (default: fasta)", default="fasta")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    return parser
def main_filter_short_seqs(args):
    # orig by rsealfon, edited by dpark
    with util.file.open_or_gzopen(args.inFile) as inf:
        with util.file.open_or_gzopen(args.outFile, 'w') as outf:
            Bio.SeqIO.write(
                [s for s in Bio.SeqIO.parse(inf, args.format)
                    if len(s) >= args.minLength],
                outf, args.output_format)
    return 0
__commands__.append(('filter_short_seqs', main_filter_short_seqs, parser_filter_short_seqs))
