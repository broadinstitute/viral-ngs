#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

# Only sometimes cache the tools/build directory
if [ -z "$TRAVIS_TAG" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $CACHE_DIR/tools_build
    ln -s $CACHE_DIR/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on tag: $TRAVIS_TAG"

fi

# Report how big things are
echo "Docker cache space usage:"
du -hs $CACHE_DIR/*
