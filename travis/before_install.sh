#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e -o pipefail

# Only sometimes cache the tools/build directory
if [ -z "$TRAVIS_TAG" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $CACHE_DIR/tools_build $CACHE_DIR/conda-tools $CACHE_DIR/conda-cache
    # ln -s $CACHE_DIR/tools_build tools/build
    #ln -s $CACHE_DIR/conda-cache tools/conda-cache
    ln -s $CACHE_DIR/conda-tools tools

    # Report how big things are
    echo "Docker cache space usage:"
    du -hs $MINICONDA_DIR $CACHE_DIR/*
else
    echo "Travis docker cache disabled for tools/build on tag: $TRAVIS_TAG"
fi
