#!/usr/bin/env python

"""Calculate the memory allocated to the process, taking account of cgroups.
Print result to stdout.
"""

import argparse
import sys
import os
import os.path
import multiprocessing

parser = argparse.ArgumentParser('Calculated memory allocated to the process')
parser.add_argument('mem_unit', choices=('mb', 'gb'), help='memory units')
parser.add_argument('mem_fraction', type=int, help='what fraction of total memory to report')
parser.add_argument('--per-cpu', dest="per_cpu", action='store_true', help='Calculate memory per-CPU.')
args = parser.parse_args()

if not (1 <= args.mem_fraction <= 100):
    raise RuntimeError("mem_fraction should be in the range [1,100]")

unit2factor = {'k': 1024, 'm': 1024*1024, 'g': 1024*1024*1024}

def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """

    cgroup_cpus = MAX_INT32
    try:
        def slurp_file(fname):
            with open(fname) as f:
                return f.read()
        def get_cpu_val(name):
            return float(slurp_file('/sys/fs/cgroup/cpu/cpu.'+name).strip())
        cfs_quota = get_cpu_val('cfs_quota_us')
        if cfs_quota > 0:
            cfs_period = get_cpu_val('cfs_period_us')
            log.debug('cfs_quota %s, cfs_period %s', cfs_quota, cfs_period)
            cgroup_cpus = max(1, int(cfs_quota / cfs_period))
    except Exception as e:
        pass

    proc_cpus = MAX_INT32
    try:
        with open('/proc/self/status') as f:
            status = f.read()
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                proc_cpus = res
    except IOError:
        pass

    log.debug('cgroup_cpus %d, proc_cpus %d, multiprocessing cpus %d',
              cgroup_cpus, proc_cpus, multiprocessing.cpu_count())
    return min(cgroup_cpus, proc_cpus, multiprocessing.cpu_count())

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

if args.per_cpu:
    mem_in_bytes = mem_in_bytes/available_cpu_count()

mem_in_units = float(mem_in_bytes) / float(unit2factor[args.mem_unit[0]])
print(int(mem_in_units * (float(args.mem_fraction) / 100.0)))
