#!/bin/bash
set -e

CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"
PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`
if [ "$PYVER" = "3" ]; then
    echo "pip installing snakemake packages (py3 only)"
    #conda install -q -y $CONDA_CHANNEL_STRING -p tools/conda-tools/default --file requirements-py3.txt python="$TRAVIS_PYTHON_VERSION"
    pip install --quiet -r requirements-py3.txt
fi

#python --version
