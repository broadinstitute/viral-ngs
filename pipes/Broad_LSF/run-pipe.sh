#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# determine the directory of this script
SCRIPT_DIRECTORY=$(dirname $0)

# load necessary Broad dotkits
source /broad/software/scripts/useuse 
reuse -q LSF
source $SCRIPT_DIRECTORY/../Broad_common/setup_dotkits.sh

# load config dirs from config.yaml. After using the conda dotkit, we should have PyYAML 
VENVDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['venv_dir']));f.close()"`
BINDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['bin_dir']));f.close()"`
DATADIR=`python -c "import yaml, os; f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['data_dir']));f.close()"`

# load Python virtual environment
source "$VENVDIR/bin/activate"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
	--jobs 100000 --immediate-submit \
        --latency-wait 20 \
	--config mode=LSF job_profiler="$BINDIR/pipes/Broad_LSF/lsf-report.py" \
	--directory . \
	--jobscript "$BINDIR/pipes/Broad_LSF/jobscript.sh" \
	--cluster $BINDIR"/pipes/Broad_LSF/cluster-submitter.py {dependencies} $DATADIR {config[log_dir]}" \
	"$@"
