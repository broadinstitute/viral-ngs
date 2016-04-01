#!/bin/bash
# Wrappers around Snakemake for use on the Broad UGER cluster

# determine the directory of this script
SCRIPT_DIRECTORY=$(dirname $(readlink --canonicalize-existing $0))

# load necessary Broad dotkits
source /broad/software/scripts/useuse 
reuse -q UGER
source $SCRIPT_DIRECTORY/../Broad_common/setup_dotkits.sh

# load config dirs from config.yaml. After using the conda dotkit, we should have PyYAML 
VENVDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['venv_dir']));f.close()"`
BINDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['bin_dir']));f.close()"`
DATADIR=`python -c "import yaml, os; f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['data_dir']));f.close()"`

# load Python virtual environment
source "$VENVDIR/bin/activate"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
    --jobs 100 --immediate-submit \
        --latency-wait 60 \
    --config mode=UGER \
    --directory . \
    --jobscript "$BINDIR/pipes/Broad_UGER/jobscript.sh" \
    --cluster $BINDIR'/pipes/Broad_UGER/cluster-submitter.py {dependencies} {config[log_dir]}' \
    "$@"
