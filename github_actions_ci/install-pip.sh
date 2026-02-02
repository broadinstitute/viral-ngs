#!/bin/bash
set -e

CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"
PYVER=`echo $GITHUB_ACTIONS_PYTHON_VERSION | cut -c 1`
if [ "$PYVER" = "3" ]; then
    echo "pip installing snakemake packages (py3 only)"
    #conda install -q -y $CONDA_CHANNEL_STRING -p tools/conda-tools/default --file requirements-py3.txt python="$GITHUB_ACTIONS_PYTHON_VERSION"
    pip install --quiet -r requirements-py3.txt
elif [ "$PYVER" = "2" ]; then
    echo "pip install py2 packages"
    #conda install -q -y $CONDA_CHANNEL_STRING -p tools/conda-tools/default --file requirements-py2.txt python="$GITHUB_ACTIONS_PYTHON_VERSION"
    pip install --quiet -r requirements-py2.txt
fi

#python --version
