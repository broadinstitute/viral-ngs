#!/usr/bin/env python3
import os
import sys
import re
import getpass
from snakemake.utils import read_job_properties

LOGDIR = sys.argv[-2]
jobscript = sys.argv[-1]
mo = re.match(r'(\S+)/snakejob\.\S+\.(\d+)\.sh', jobscript)
assert mo
sm_tmpdir, sm_jobid = mo.groups()
props = read_job_properties(jobscript)

# Blacklist problematic nodes; this list is stored as filenames
# in /broad/hptmp/[username]/blacklisted-nodes/
whoami = getpass.getuser()
blacklisted_node_dir = os.path.join("/broad/hptmp", whoami, "blacklisted-nodes")
if not os.path.exists(blacklisted_node_dir):
    os.makedirs(blacklisted_node_dir)
def hard_blacklist_node(node):
    blacklist_path = os.path.join(blacklisted_node_dir, node)
    with open(blacklist_path, 'a'):
        os.utime(blacklist_path, None)
# Always blacklist 'sgi1'; it cannot perform basic operations like
# allocating memory
hard_blacklist_node('sgi1')
blacklisted_nodes = os.listdir(blacklisted_node_dir)

# set up job name, project name
jobname = "{rule}-{jobid}".format(rule=props["rule"], jobid=sm_jobid)
if props["params"].get("logid"):
    jobname = "{rule}-{id}".format(rule=props["rule"], id=props["params"]["logid"])
cmdline = "qsub -N {jobname} -cwd -r y ".format(jobname=jobname)

# log file output
cmdline += "-o {logdir} -j y ".format(logdir=LOGDIR)

# pass memory resource request to cluster
mem = props.get('resources', {}).get('mem_mb')
threads = props.get('resources', {}).get('threads')

if mem:
    mem = int(mem)
    threads = int(threads) or 1 # only used here for the calculation of memory per-core
    # on UGER, memory requests are per-core (according to BITS as of Sept. 6, 2016)
    mem_per_core = round(float(mem)/float(threads), 2)
    # increase memory requested to reflect new JVM-UGER changes requiring greater memory headroom
    mem_per_core = round(mem_per_core*1.1,2)
    if blacklisted_nodes:
        # Pass h= as the hostname parameter; it accepts a regex, so
        # invert the match to blacklist hostnames
        cmdline += " -l h_vmem={}M,h_rss={}M,h='!({})' ".format(
            mem_per_core, round(1.2 * mem_per_core, 2),
            '|'.join(blacklisted_nodes))
    else:
        cmdline += " -l h_vmem={}M,h_rss={}M ".format(
            mem_per_core, round(1.2 * mem_per_core, 2))
    if mem >= 15 or (threads and threads >= 4):
        cmdline += ' -R y '
elif blacklisted_nodes:
    # Pass h= as the hostname parameter; it accepts a regex, so
    # invert the match to blacklist hostnames
    cmdline += " -l h='!({})' ".format('|'.join(blacklisted_nodes))

if threads:
    cmdline += ' -pe smp {} -binding linear:{} '.format(int(threads), int(threads))

# rule-specific UGER parameters (e.g. queue)
cmdline += props["params"].get("UGER", "") + " "

# figure out job dependencies
dependencies = set(sys.argv[1:-2])
if dependencies:
    cmdline += "-hold_jid '{}' ".format(",".join(dependencies))

# the actual job
cmdline += jobscript

# the success file
cmdline += " %s/%s.jobfinished" % (sm_tmpdir, sm_jobid)

# the part that strips bsub's output to just the job id
cmdline += r" | tail -1 | cut -f 3 -d \ "

# call the command
os.system(cmdline)
