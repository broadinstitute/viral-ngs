#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# load config dirs from config.yaml
VENVDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["venv_dir"]));f.close()'`
BINDIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["bin_dir"]));f.close()'`

source "$BINDIR/pipes/Broad_UGER/setup_dotkits.sh"

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
