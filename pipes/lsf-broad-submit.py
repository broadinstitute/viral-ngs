#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
props = read_job_properties(jobscript)

cmdline = "bsub -N -P {proj_name} -J {job_name} ".format(
    proj_name='viral_ngs',
    job_name=props["rule"])

cmdline += props["params"].get("LSF","") + " "

# figure out job dependencies
dependencies = sys.argv[1:-1]
if dependencies:
    cmdline += "-w '{}' ".format(" && ".join(dependencies))

# the actual job
cmdline += jobscript

# the part that strips bsub's output to just the job id
cmdline += " | cut -f 2 -d \< | cut -f 1 -d \>"

# call the command
os.system(cmdline)
