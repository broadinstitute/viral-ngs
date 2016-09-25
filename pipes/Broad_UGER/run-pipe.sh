#!/bin/bash
# Wrappers around Snakemake for use on the Broad UGER cluster

# determine the directory of this script
SCRIPT_DIRECTORY=$(dirname $(readlink --canonicalize-existing $0))

# load necessary Broad dotkits
source /broad/software/scripts/useuse
reuse -q UGER

reuse -q Python-3.4
# load config dirs from config.yaml. After using the conda dotkit, we should have PyYAML 
CONDAENVDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['conda_env_dir']));f.close()"`
MINICONDADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["miniconda_dir"]));f.close()'`
BINDIR=`python -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['bin_dir']));f.close()"`
DATADIR=`python -c "import yaml, os; f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['data_dir']));f.close()"`
unuse  Python-3.4

export PATH="$MINICONDADIR/bin:$PATH"

# load conda environment
source activate "$CONDAENVDIR"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
    --jobs 90 --immediate-submit \
    --force-use-threads \
    --latency-wait 60 \
    --config mode=UGER \
    --directory . \
    --jobscript "$BINDIR/pipes/Broad_UGER/jobscript.sh" \
    --cluster $BINDIR'/pipes/Broad_UGER/cluster-submitter.py {dependencies} {config[log_dir]}' \
    "$@"
