#!/bin/bash
set -e -o pipefail

CONDA_ENV="$(pwd)/tools/conda-tools/default"

# Set to conda's java
export JAVA_HOME="$CONDA_ENV/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`

if [ ! -d $CONDA_ENV ]; then
	conda create -y -m -p $CONDA_ENV python="$TRAVIS_PYTHON_VERSION"
fi

conda install -y --override-channels \
	-c broad-viral -c r -c bioconda -c conda-forge -c defaults \
	--file requirements-conda.txt \
	--file requirements-conda-tests.txt \
	--file requirements-py$PYVER.txt \
	-p $CONDA_ENV

conda list -p $CONDA_ENV

conda clean --all --yes # clean temp/cache files to reduce Travis cache size
