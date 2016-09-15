#!/usr/bin/env python
''' This script contains utilities for analyzing LSF reports for
    data on runtime and performance.
'''

__author__ = "dpark@broadinstitute.org"

import argparse
import re
import time
import os
import os.path
import sys


def read_lsf_logfile(infname):
    out = {'logfile': infname}
    num_dash_lines = 0
    with open(infname, 'rU', encoding='latin-1') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('Subject:'):
                mo = re.match(r'Subject: Job (\d+)', line)
                out['job_id'] = mo.group(1)
            elif line.startswith('Job <'):
                mo = re.match(r'Job <(\S+)> was submitted from host', line)
                out['job_name'] = mo.group(1)
                if '-' in out['job_name']:
                    out['job_prefix'], out['job_suffix'] = out['job_name'].split('-', 1)
            elif line.startswith('Job was executed on host'):
                mo = re.match(r'Job was executed on host\(s\) <(\S+?)>, in queue <(\w+)>', line)
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
                        if line.startswith('TERM'):
                            line = line.split(':')[0]
                        out['status'] = line
                    elif ':' in line:
                        k, v = [s.strip() for s in line.split(':')]
                        out[k] = v
    if 'start_time' in out and 'end_time' in out:
        out['run_time'] = time.mktime(time.strptime(out['end_time'])) \
            - time.mktime(time.strptime(out['start_time']))
    return out


def read_all_logfiles(dirname):
    header = [
        'job_id', 'job_name', 'job_prefix', 'job_suffix', 'queue', 'exec_host', 'status', 'run_time', 'start_time',
        'end_time', 'CPU time', 'Max Memory', 'Max Swap', 'Max Processes', 'Max Threads', 'logfile'
    ]
    yield header
    for fname in os.listdir(dirname):
        try:
            row = read_lsf_logfile(os.path.join(dirname, fname))
        except:
            print("Error parsing " + fname)
            raise
        yield [str(row.get(h, '')) for h in header]


def parser_report():
    parser = argparse.ArgumentParser(description="Read a directory full of LSF log files and produce a tabular report.")
    parser.add_argument("log_dir", help="Input directory of LSF log files")
    parser.add_argument("outFile", help="Output report file")
    return parser


def main_report(args):
    with open(args.outFile, 'wt') as outf:
        for row in read_all_logfiles(args.log_dir):
            outf.write('\t'.join(row) + '\n')
    return 0


if __name__ == '__main__':
    argv = sys.argv[1:]
    report_parser = parser_report()
    if len(argv) == 0:
        report_parser.print_help()
    else:
        report_args = report_parser.parse_args(argv)
        main_report(report_args)
