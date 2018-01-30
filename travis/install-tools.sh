#!/bin/bash
set -e -o pipefail

# Set to conda's java
export JAVA_HOME="$(pwd)/tools/conda-tools/default/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`

conda create -y -m -p tools/conda-tools/default python="$TRAVIS_PYTHON_VERSION"

conda install -y --override-channels \
	-c broad-viral -c r -c bioconda -c conda-forge -c defaults \
	--file requirements-conda.txt \
	--file requirements-conda-tests.txt \
	--file requirements-py$PYVER.txt \
	-p tools/conda-tools/default

conda list -p tools/conda-tools/default

conda clean --all --yes # clean temp/cache files to reduce Travis cache size
