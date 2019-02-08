#!/usr/bin/env python3
import os
import sys
import re
from snakemake.utils import read_job_properties

LOGDIR = sys.argv[-2]
DATADIR = sys.argv[-3]
jobscript = sys.argv[-1]
mo = re.match(r'(\S+)/snakejob\.\S+\.(\d+)\.sh', jobscript)
assert mo
sm_tmpdir, sm_jobid = mo.groups()
props = read_job_properties(jobscript)

# set up job name, project name
jobname = "{rule}-{jobid}".format(rule=props["rule"], jobid=sm_jobid)
if props["params"].get("logid"):
    jobname = "{rule}-{id}".format(rule=props["rule"], id=props["params"]["logid"])

# -E is a pre-exec command, that reschedules the job if the command fails
#   in this case, if the data dir is unavailable (as may be the case for a hot-mounted file path)
cmdline = 'bsub -P {proj_name} -J {jobname} -r -E "ls {datadir}" '.format(proj_name='viral_ngs', jobname=jobname, datadir=DATADIR)

# log file output
if "-N" not in props["params"].get("LSF", ""):
    cmdline += "-oo {logdir}/LSF-{jobname}.txt ".format(logdir=LOGDIR, jobname=jobname)

# pass memory resource request to LSF
mem = int(props.get('resources', {}).get('mem_mb',1000))/1000
if mem:
    cmdline += '-R "rusage[mem={}]" -M {} '.format(mem, 2 * int(mem))

# rule-specific LSF parameters (e.g. queue, runtime)
cmdline += props["params"].get("LSF", "") + " "

# figure out job dependencies
dependencies = set(sys.argv[1:-3])
if dependencies:
    cmdline += "-w '{}' ".format(" && ".join(dependencies))

# the actual job
cmdline += jobscript

# the success file
cmdline += " %s/%s.jobfinished" % (sm_tmpdir, sm_jobid)

# the part that strips bsub's output to just the job id
cmdline += r" | tail -1 | cut -f 2 -d \< | cut -f 1 -d \>"

# call the command
os.system(cmdline)
