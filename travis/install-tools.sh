#!/bin/bash
set -e -o pipefail

CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"
CONDA_ENV="$(pwd)/tools/conda-tools/default"

# Set to conda's java
export JAVA_HOME="$CONDA_ENV/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`

if [ ! -d $CONDA_ENV ]; then
	conda create -y -m --quiet -p $CONDA_ENV python="$TRAVIS_PYTHON_VERSION"
fi

conda install -y --quiet \
	$CONDA_CHANNEL_STRING \
	--file requirements-conda.txt \
	--file requirements-conda-tests.txt \
	--file requirements-py$PYVER.txt \
	-p $CONDA_ENV

conda clean --all --yes # clean temp/cache files to reduce Travis cache size

conda list -p $CONDA_ENV
