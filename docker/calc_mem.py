#!/usr/bin/env python

"""Calculate the memory allocated to the process, taking account of cgroups.
Print result to stdout.
"""

import argparse
import sys
import os
import os.path

parser = argparse.ArgumentParser('Calculated memory allocated to the process')
parser.add_argument('mem_unit', choices=('mb', 'gb'), help='memory units')
parser.add_argument('mem_fraction', type=int, help='what fraction of total memory to report')
args = parser.parse_args()

if not (1 <= args.mem_fraction <= 100):
    raise RuntimeError("mem_fraction should be in the range [1,100]")

unit2factor = {'k': 1024, 'm': 1024*1024, 'g': 1024*1024*1024}

def mem_from_proc_meminfo():
    """Return the total memory, in bytes, as given by /proc/meminfo"""
    with open('/proc/meminfo') as f:
        for line in f:
            if line.startswith('MemTotal:'):
                parts = line.strip().split()
                val, unit = parts[1:3]
                unit_factor = unit2factor[unit[0].lower()]
                return int(val) * unit_factor
    raise RuntimeError('Could not get MemTotal from /proc/meminfo')

def mem_from_cgroups():
    """Return the total memory, in bytes, as given by cgroups (or sys.maxsize if not given)"""
    cgroups_memlimit_fname = '/sys/fs/cgroup/memory/memory.limit_in_bytes'
    if os.path.isfile(cgroups_memlimit_fname):
        with open(cgroups_memlimit_fname) as f:
            val = f.read().strip()
            return int(val) * unit2factor.get(val[-1], 1)

    return sys.maxsize

mem_in_bytes = min(mem_from_proc_meminfo(), mem_from_cgroups())
mem_in_units = float(mem_in_bytes) / float(unit2factor[args.mem_unit[0]])
print(int(mem_in_units * (float(args.mem_fraction) / 100.0)))
