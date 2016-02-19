#!/bin/bash
# Wrappers around Snakemake for use on the Broad UGER cluster

# determine the directory of this script
SCRIPT_DIRECTORY=$(dirname $0)

# load necessary Broad dotkits
source /broad/software/scripts/useuse 
reuse -q UGER
source $SCRIPT_DIRECTORY/../Broad_common/setup_dotkits.sh

# uses the first argument as the config file path, is specified
if [[ ! -z "$1" && "$1" != " " ]]; then
    CONFIG_FILE=$1
else
    # otherwise the config file is assumed to be "config.yaml" in the cwd
    CONFIG_FILE="config.yaml"
fi

# resolve the config file path in full
CONFIG_FILE=`python -c "import os; print( os.path.realpath(os.path.expanduser(\"$CONFIG_FILE\"))) "`

# if the config file does not exist
# it is either not in the cwd, or what the user passed in does not exist
if [[ ! -f $CONFIG_FILE ]]; then
    echo "Config file does not exist: $CONFIG_FILE"
    echo "   Usage: $(basename $0) [path/to/config.yaml]"
    echo "   A file called 'config.yaml' must exist in the current directory, or be passed in."
    exit 1
fi

# load config dirs from config.yaml. After using the conda dotkit, we should have PyYAML 
VENVDIR=`python -c "import yaml, os;f=open(\"$CONFIG_FILE\");print(os.path.realpath(yaml.safe_load(f)['venv_dir']));f.close()"`
BINDIR=`python -c "import yaml, os;f=open(\"$CONFIG_FILE\");print(os.path.realpath(yaml.safe_load(f)['bin_dir']));f.close()"`
DATADIR=`python -c "import yaml, os; f=open(\"$CONFIG_FILE\");print(os.path.realpath(yaml.safe_load(f)['data_dir']));f.close()"`

# load Python virtual environment
source "$VENVDIR/bin/activate"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
    --jobs 100 --immediate-submit \
        --latency-wait 20 \
    --config mode=UGER \
    --directory . \
    --jobscript "$BINDIR/pipes/Broad_UGER/jobscript.sh" \
    --cluster $BINDIR'/pipes/Broad_UGER/cluster-submitter.py {dependencies} {config[log_dir]}' \
    "$@"
