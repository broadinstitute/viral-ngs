#!/bin/bash
set -e -o pipefail

# Set to conda's java
export JAVA_HOME="$(pwd)/tools/conda-tools/default/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`

for i in $(seq 3); do
  conda create --quiet -y -m -c broad-viral -c r -c bioconda -c conda-forge -c defaults -p tools/conda-tools/default --file requirements-conda.txt --file requirements-conda-tests.txt --file requirements-py$PYVER.txt python="$TRAVIS_PYTHON_VERSION" && break
  sleep 5
done

conda list

conda clean --all --yes # clean temp/cache files to reduce Travis cache size
