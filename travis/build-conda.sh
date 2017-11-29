#!/bin/bash
# This script performs various packing and deployment operations.

set -e -o pipefail

# === Build conda package
CONDA_PACKAGE_OUTDIR=packaging/conda-packages
echo "Rendering and building conda package..."
echo "Python binary: $(which python)"
echo "Python version: $(python --version)"
# Render recipe from template and dependency files, setting the tag as the current version
# if this is a tag build+upload, otherwise just test building

# if the ANACONDA_TOKEN is defined (not on an external branch)
if [ -z "$ANACONDA_TOKEN" ]; then
    echo "ANACONDA_TOKEN is not defined. Conda package upload is only supported for branches on the original repository."
    exit 1
fi

conda config --set anaconda_upload yes
if [ -n "$TRAVIS_TAG" ]; then
    # This is an official release

    # render and build the conda package
    python packaging/conda-recipe/render-recipe.py "$TRAVIS_TAG" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt
    CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" packaging/conda-recipe/viral-ngs

else
    # This is a development build

    # make a directory to hold the built conda package
    mkdir -p CONDA_PACKAGE_OUTDIR
    
    if [ -z "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
        # if a commit is being pushed, TRAVIS_PULL_REQUEST_BRANCH is empty
        BRANCH_NAME="$TRAVIS_BRANCH"
    else
        BRANCH_NAME="$TRAVIS_PULL_REQUEST_BRANCH"
    fi

    SANITIZED_BRANCH_NAME="$(echo $BRANCH_NAME | sed 's/-/_/g')"
    CONDA_PKG_VERSION="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1+dev-/' | sed 's/-/_/g')_$(echo $SANITIZED_BRANCH_NAME)"
    echo "Building conda package version $CONDA_PKG_VERSION"

    # render and build the conda package
    python packaging/conda-recipe/render-recipe.py "$CONDA_PKG_VERSION" --package-name "viral-ngs-dev" --download-filename "$TRAVIS_COMMIT" --build-reqs requirements-conda.txt --run-reqs requirements-conda.txt --py3-run-reqs requirements-py3.txt --py2-run-reqs requirements-py2.txt --test-reqs requirements-conda-tests.txt
    CONDA_PERL=5.22.0 conda build -c broad-viral -c r -c bioconda -c conda-forge -c defaults --python "$TRAVIS_PYTHON_VERSION" --token "$ANACONDA_TOKEN" --output-folder "$CONDA_PACKAGE_OUTDIR" packaging/conda-recipe/viral-ngs

fi
