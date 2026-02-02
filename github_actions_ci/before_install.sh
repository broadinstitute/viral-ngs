#!/bin/bash
# This script primarily enables or disables the CI dependency
# cache depending on whether we're on the master branch or not.
set -e -o pipefail

# Only sometimes cache the tools/build directory
if [ -z "$GITHUB_ACTIONS_TAG" ]; then
    echo "CI docker caches allowed for branch $GITHUB_ACTIONS_BRANCH"
    rm -rf tools/build
    mkdir -p $CACHE_DIR/tools_build $CACHE_DIR/conda-tools $CACHE_DIR/conda-cache
    # ln -s $CACHE_DIR/tools_build tools/build
    #ln -s $CACHE_DIR/conda-cache tools/conda-cache
    ln -s $CACHE_DIR/conda-tools tools

    # Report how big things are
    echo "Docker cache space usage:"
    du -hs $MINICONDA_DIR $CACHE_DIR/*
else
    echo "CI docker cache disabled for tools/build on tag: $GITHUB_ACTIONS_TAG"
fi
