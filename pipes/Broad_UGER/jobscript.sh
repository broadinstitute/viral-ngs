#!/bin/sh
# properties = {properties}
# this is identical to the default jobscript with the exception of the exit code

source /broad/software/scripts/useuse
use UGER
use .anaconda3-5.0.1
CONDAENVDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["conda_env_dir"]));f.close()'`
MINICONDADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["miniconda_dir"]));f.close()'`
BINDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["bin_dir"]));f.close()'`
DATADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["data_dir"]));f.close()'`
unuse .anaconda3-5.0.1

export PATH="$MINICONDADIR/bin:$PATH"

# load Python virtual environment
source activate "$CONDAENVDIR"


# As a simple solution to transient UGER problems, maintain a list of blacklisted
# nodes. When a node fails a check, add its hostname to the list and exit 99
# to reschedule the job. Blacklist these nodes via qsub
# Maintain the list of blacklisted nodes as filenames in /broad/hptmp/$(whoami)/blacklisted-nodes
BLACKLISTED_NODES="/broad/hptmp/$(whoami)/blacklisted-nodes"
mkdir -p $BLACKLISTED_NODES

# Cleanup blacklisted nodes if they have not been touched in a day
find $BLACKLISTED_NODES -name "*" -type f -mmin +2880 -delete

# Specify whether to require that a node have the NFS share mounted
REQUIRE_NFS_SHARE_MOUNTED=true

# Specify whether to require that a node have >1G of shared memory
# (/dev/shm) available (e.g., snakemake requires a little for its
# use of multiprocessing.Lock() to setup logging)
REQUIRE_AVAILABLE_SHARED_MEMORY=true

# Specify whether to require that a node can run java without
# error (e.g., for some time the node 'sgi1' could not run java)
REQUIRE_JAVA_TO_RUN=true

# Perform checks on node and decide whether to blacklist
if [[ "$REQUIRE_NFS_SHARE_MOUNTED" = true ]] && ! $(ls "$DATADIR" &> /dev/null); then
    # Listing the data directory fails since the node does not have the
    # NFS share mounted
    echo "Host '$(hostname)' does not have NFS share mounted. Retrying.." 1>&2
    touch "$BLACKLISTED_NODES/$(hostname)"
    exit 99
fi
if [[ "$REQUIRE_AVAILABLE_SHARED_MEMORY" = true ]] && [[ $(df -k /dev/shm | tail -n 1 | awk '{{print $4}}') -lt 1000000 ]]; then
    # There is too little shared memory available
    echo "Host '$(hostname)' has full or near-full /dev/shm. Retrying.." 1>&2
    touch "$BLACKLISTED_NODES/$(hostname)"
    exit 99
fi
if [[ "$REQUIRE_JAVA_TO_RUN" = true ]] && ! java -version; then
    # Java does not run successfully
    echo "Host '$(hostname)' cannot run java. Retrying.." 1>&2
    touch "$BLACKLISTED_NODES/$(hostname)"
    exit 99
fi


echo -e "JOB ID\t$JOB_ID"
echo "=============================="

{exec_job}
EXIT_STATUS=$?

# Report resource consumption because it's not reported by default
echo "------------------------------"
qstat -j $JOB_ID | grep '^usage'

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null
if [[ $? -eq 0 ]]; then
    exit 0
else
    if [[ "{workflow.immediate_submit}" -eq "True" ]]; then
        exit 100
    else
        exit $EXIT_STATUS
    fi
fi
