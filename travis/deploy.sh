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

# if the ANACONDA_TOKEN is defined (not on an external branch)
if [ ! -z "$ANACONDA_TOKEN" ]; then
    echo "Running $SCRIPTPATH/package-conda.sh"
    # Render recipe from template and dependency files, setting the tag as the current version
    packaging/conda-recipe/render-recipe.py "$PKG_VERSION" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt && \
        conda build --python "$TRAVIS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" packaging/conda-recipe/viral-ngs && \
        ./travis/trigger-tests-in-other-repo.sh
        # check the exit code of conda build, and if successful,
        # trigger the viral-ngs-deploy repository to test/build the docker container
else
    echo "ANACONDA_TOKEN is not defined. Conda package upload is only supported for branches on the original repository."
fi
