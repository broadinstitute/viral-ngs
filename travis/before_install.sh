#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

# Only sometimes cache the tools/build directory
if [ -z "$TRAVIS_TAG" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $CACHE_DIR/tools_build $CACHE_DIR/conda-tools $CACHE_DIR/conda-cache
    ln -s $CACHE_DIR/tools_build tools/build
    #ln -s $CACHE_DIR/conda-cache tools
    #ln -s $CACHE_DIR/conda-tools tools

    # Report how big things are
    echo "Docker cache space usage:"
    du -hs $MINICONDA_DIR $PIP_DIR $CACHE_DIR/*
else
    echo "Travis docker cache disabled for tools/build on tag: $TRAVIS_TAG"
fi

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # used to correct this error on Travis:
    #   /Users/travis/build.sh: line 159: shell_session_update: command not found
    # Under OSX, some versions of ruby seem to cause this error.
    # See: https://github.com/travis-ci/travis-ci/issues/6307
    rvm get stable
fi
