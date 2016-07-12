#!/bin/sh
# properties = {properties}
# this is identical to the default jobscript with the exception of the exit code

source /broad/software/scripts/useuse
use Python-3.4
CONDAENVDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["conda_env_dir"]));f.close()'`
MINICONDADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["miniconda_dir"]));f.close()'`
BINDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["bin_dir"]));f.close()'`
DATADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["data_dir"]));f.close()'`
unuse Python-3.4

export PATH="$MINICONDADIR/bin:$PATH"

# load Python virtual environment
source activate "$CONDAENVDIR/bin"

# if listing the data directory fails, exit 99 to reschedule the job
# since the node with the job doesn't have the NFS share mounted.
# Hopefully it will be sent to a node that has it mounted.
# For this to work, re-run must be set to y via "qsub -r y"
ls "$DATADIR" || exit 99

echo $JOB_ID
echo "=============================="

{exec_job}

# Report resource consumption because it's not reported by default
echo "------------------------------"
qstat -j $JOB_ID | grep '^usage'

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null && exit 0 || exit 100
