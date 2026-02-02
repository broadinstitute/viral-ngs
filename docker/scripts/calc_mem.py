#!/usr/bin/env python

"""Calculate the memory allocated to the process, taking account of cgroups.
Print result to stdout.
"""

import math
import argparse
import sys
import os
import re
import logging
import os.path
import multiprocessing

import psutil

#from util.misc import available_cpu_count # use the version of available_cpu_count() from viral-core/util/misc.py

log = logging.getLogger(__name__)

parser = argparse.ArgumentParser('Calculated memory allocated to the process')
parser.add_argument('mem_unit', choices=('b', 'kb', 'mb', 'gb'), help='memory units')
parser.add_argument('mem_fraction', type=int, help='what fraction of total memory to report')
parser.add_argument('--per-cpu', dest="per_cpu", action='store_true', help='Calculate memory per-CPU.')
args = parser.parse_args()

if not (1 <= args.mem_fraction <= 100):
    raise RuntimeError("mem_fraction should be in the range [1,100]")

unit2factor = {'b': 1, 'k': 1024, 'm': 1024*1024, 'g': 1024*1024*1024}
MAX_INT32 = (2 ** 31)-1

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
        def _load(path, encoding="utf-8"):
            """ Loads a file content """
            with open(path, 'r', encoding=encoding, newline="") as handle:
                tmp = handle.read()
            return tmp

        # cgroup CPU count determination (w/ v2) adapted from:
        #   https://github.com/conan-io/conan/blob/2.9.2/conan/tools/build/cpu.py#L31-L54
        #
        # see also:
        #   https://docs.kernel.org/scheduler/sched-bwc.html

        # This is necessary to determine docker cpu_count
        cfs_quota_us = cfs_period_us = 0
        # cgroup v2
        if os.path.exists("/sys/fs/cgroup/cgroup.controllers"):
            cpu_max = _load("/sys/fs/cgroup/cpu.max").split()
            if cpu_max[0] != "max":
                if len(cpu_max) == 1:
                    cfs_quota_us, cfs_period_us = int(cpu_max[0]), 100_000
                else:
                    cfs_quota_us, cfs_period_us = map(int, cpu_max)
        # cgroup v1
        else:
            cfs_quota_us = int(_load("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"))
            cfs_period_us = int(_load("/sys/fs/cgroup/cpu/cpu.cfs_period_us"))

        log.debug('cfs_quota_us %s, cfs_period_us %s', cfs_quota_us, cfs_period_us)
        if cfs_quota_us > 0 and cfs_period_us > 0:
            cgroup_cpus = max(1, int(math.ceil(cfs_quota_us / cfs_period_us)))
    except Exception as e:
        pass

    proc_cpus = MAX_INT32
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', _load('/proc/self/status'))
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
    # list of potential cgroup paths to max mem info
    #   see:
    #     (cgroup v1) https://www.kernel.org/doc/Documentation/cgroup-v1/memory.txt
    #     (cgroup v2) https://www.kernel.org/doc/Documentation/cgroup-v2.txt
    cgroups_memlimit_fnames = [
                                '/sys/fs/cgroup/memory/memory.limit_in_bytes', # cgroup v1
                                '/sys/fs/cgroup/memory.max' # cgroup v2
                             ]
    # try the various potential cgroup memory info paths
    for cgroups_memlimit_fname in cgroups_memlimit_fnames:
        if os.path.isfile(cgroups_memlimit_fname):
            with open(cgroups_memlimit_fname) as f:
                val = f.read().strip()
                if val != "max":
                    return int(val) * unit2factor.get(val[-1], 1)

    return sys.maxsize

def mem_from_psutil(metric_name="total"):
    """ Use psutil to get a memory metric by name in a cross-platform way
        Returning sys.maxsize (obviously wrong large value) 
        in the event the value cannot be obtained.

        For available metrics, see:
          https://psutil.readthedocs.io/en/latest/#psutil.virtual_memory
    """
    mem_info = psutil.virtual_memory()

    return int(getattr(mem_info,metric_name,sys.maxsize))

# of the memory values obtained, use the smallest value
# this results in obviously-wrong values obtained from sys.maxsize
# in mem_from_cgroups() or mem_from_psutil() falling in precedence 
mem_in_bytes = min(
                    mem_from_psutil(),
                    mem_from_proc_meminfo(),
                    mem_from_cgroups()
               )

if args.per_cpu:
    mem_in_bytes = mem_in_bytes/available_cpu_count()

mem_in_units = float(mem_in_bytes) / float(unit2factor[args.mem_unit[0]])
print(int(mem_in_units * (float(args.mem_fraction) / 100.0)))
