#!/bin/bash
# This script performs various packing and deployment operations.

set -e -o pipefail

# === Build conda package
CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"
CONDA_PACKAGE_OUTDIR=packaging/conda-packages
echo "Rendering and building conda package..."
echo "Python binary: $(which python)"
echo "Python version: $(python --version 2>&1)"
# Render recipe from template and dependency files, setting the tag as the current version
# if this is a tag build+upload, otherwise just test building

# if the ANACONDA_TOKEN is defined (not on an external branch)
if [ -z "$ANACONDA_TOKEN" ]; then
    echo "ANACONDA_TOKEN is not defined. Conda package upload is only supported for branches on the original repository."
    exit 1
fi

conda config --set anaconda_upload no
if [ -n "$GITHUB_ACTIONS_TAG" ]; then
    # This is an official release, upload it
    conda config --set anaconda_upload yes

    # render and build the conda package
    echo "Rendering recipe..."
    python packaging/conda-recipe/render-recipe.py "$GITHUB_ACTIONS_TAG" --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt # --build-reqs requirements-conda.txt
    echo "Building recipe..."
    CONDA_PERL=5.26 conda build $CONDA_CHANNEL_STRING --python "$GITHUB_ACTIONS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" packaging/conda-recipe/viral-ngs

else
    # This is a development build

    # make a directory to hold the built conda package
    mkdir -p CONDA_PACKAGE_OUTDIR
    
    if [ -z "$GITHUB_ACTIONS_PULL_REQUEST_BRANCH" ]; then
        # if a commit is being pushed, GITHUB_ACTIONS_PULL_REQUEST_BRANCH is empty
        BRANCH_NAME="$GITHUB_ACTIONS_BRANCH"
    else
        BRANCH_NAME="$GITHUB_ACTIONS_PULL_REQUEST_BRANCH"
    fi

    SANITIZED_BRANCH_NAME="$(echo $BRANCH_NAME | sed 's/-/_/g')"
    CONDA_PKG_VERSION="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-?/$1+dev/' | sed 's/-/_/g')_$(echo $SANITIZED_BRANCH_NAME)"
    echo "Building conda package version $CONDA_PKG_VERSION"

    # render and build the conda package
    echo "Rendering recipe..."
    python packaging/conda-recipe/render-recipe.py "$CONDA_PKG_VERSION" --package-name "viral-ngs-dev" --download-filename "$GITHUB_ACTIONS_COMMIT" --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt --build-reqs requirements-conda.txt
    echo "Building recipe..."
    CONDA_PERL=5.26 conda build $CONDA_CHANNEL_STRING --python "$GITHUB_ACTIONS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" --output-folder "$CONDA_PACKAGE_OUTDIR" packaging/conda-recipe/viral-ngs
fi
