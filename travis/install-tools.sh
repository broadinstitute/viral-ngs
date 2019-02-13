#!/bin/bash
set -e -o pipefail #exit on error or within pipes

CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"
CONDA_ENV="$(pwd)/tools/conda-tools/default"

# Set to conda's java
export JAVA_HOME="$CONDA_ENV/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`

if [ ! -d $CONDA_ENV ]; then
    conda create -y -m --quiet -p $CONDA_ENV python="$TRAVIS_PYTHON_VERSION" || exit $?
fi

# commented out custom retry logic for now
# until we establish the built-in conda retry logic
# works as intended when we retry more than the default
# 3 times
# remote_max_retries is set in install-conda.sh
#set +e # disable exit on error
#RETRIES=0
#RETRY_LIMIT=3
#until #conda...
conda install -y --quiet \
    $CONDA_CHANNEL_STRING \
    --file requirements-conda.txt \
    --file requirements-conda-tests.txt \
    --file requirements-py$PYVER.txt \
    -p $CONDA_ENV #; do
#        let RETRIES++
#        if [ "$RETRIES" -gt "$RETRY_LIMIT" ]; then
#            break
#        fi
#        echo "Conda install failed, likely due to transient server errors. Retrying (${RETRIES}/${RETRY_LIMIT})..."
#done
#set -e #(re)enable exit on error

conda clean --all --yes # clean temp/cache files to reduce Travis cache size

conda list -p $CONDA_ENV
