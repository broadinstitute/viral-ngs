#!/bin/bash
# Wrappers around Snakemake for execution on a single instance

# determine the directory of this script
SCRIPT_DIRECTORY=$(dirname $(readlink --canonicalize-existing $0))

# if a conda environment is active, deactivate it
# if [[ ! -z "${CONDA_PREFIX}" ]]; then
#     echo "deactivating env: $CONDA_PREFIX"
#     source deactivate
# fi

python_check=$(hash python &> /dev/null || hash python3 &> /dev/null)
if [ $? -ne 0 ]; then
    echo "It looks like Python is not installed. Exiting."
    if [[ $sourced -eq 0 ]]; then
        exit 1
    else
        return 1
    fi
fi

python3_check=$(hash python3 &> /dev/null)
if [ $? -eq 0 ]; then
    python_to_use="$(which python3)"
fi

$python_to_use --version

# load config dirs from config.yaml. After using the conda dotkit, we should have PyYAML 
CONDAENVDIR=`$python_to_use -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['conda_env_dir']));f.close()"`
MINICONDADIR=`$python_to_use -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["miniconda_dir"]));f.close()'`
BINDIR=`$python_to_use -c "import yaml, os;f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['bin_dir']));f.close()"`
DATADIR=`$python_to_use -c "import yaml, os; f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['data_dir']));f.close()"`
LOGDIR=`$python_to_use -c "import yaml, os; f=open(\"config.yaml\");print(os.path.realpath(yaml.safe_load(f)['log_dir']));f.close()"`

#export PATH="$MINICONDADIR/bin:$PATH"

# load conda environment
#source activate "$CONDAENVDIR"

ARGS=""
# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --rerun-incomplete --keep-going --nolock \
          $ARGS \
          --latency-wait 20 \
          --directory . \
          --resources mem_mb=$(expr $(cat /proc/meminfo | grep MemTotal | awk '{print $2}') / 1000) \
          --cores $(expr $(grep -c ^processor /proc/cpuinfo) - 1) \
          "$@" | tee "$LOGDIR/snakemake_$(date +%F_%s).log"