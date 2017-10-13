#!/bin/bash

# This script performs various packing and deployment operations.
# It assumes it will be caused as a deploy hook of TravisCI; ex.:
#
# deploy:
#   provider: script
#   script: travis/deploy.sh $TRAVIS_TAG
#   on:
#     tags: true
#     all_branches: master


# way to get the absolute path to this script that should
# work regardless of whether or not this script has been sourced
# Find original directory of bash script, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        SOURCE="$(readlink "$SOURCE")"
    else
        SOURCE="$(readlink -f "$SOURCE")"
    fi
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT=$SOURCE
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$SCRIPT_DIRNAME" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"

PKG_VERSION=$1

# === Build conda package and upload

echo "Python binary: $(which python)"
echo "Python version: $(python --version)"

# If this is a PR, on the master branch, or is a tag, render and build the conda package. If it is a tag, also upload to anaconda.org
if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then
    echo "Rendering and building conda package..."
    # Render recipe from template and dependency files, setting the tag as the current version
    # if this is a tag build+upload, otherwise just test building
    if [ -n "$TRAVIS_TAG" ]; then
         # if the ANACONDA_TOKEN is defined (not on an external branch)
        if [ ! -z "$ANACONDA_TOKEN" ]; then
            conda config --set anaconda_upload yes
            python packaging/conda-recipe/render-recipe.py "$PKG_VERSION" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt && \
            CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" packaging/conda-recipe/viral-ngs && \
            ./travis/trigger-tests-in-other-repo.sh
            # check the exit code of conda build, and if successful,
            # trigger the viral-ngs-deploy repository to test/build the docker container
        else
            echo "ANACONDA_TOKEN is not defined. Conda package upload is only supported for branches on the original repository."
        fi
    else
        python packaging/conda-recipe/render-recipe.py "0.0.0" --download-filename "$TRAVIS_BRANCH" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt && \
        CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --no-anaconda-upload packaging/conda-recipe/viral-ngs
    fi
fi  

