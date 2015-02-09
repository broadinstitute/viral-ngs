#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

# Only sometimes cache the tools/build directory
if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $HOME/virtualenv/tools_build
    ln -s $HOME/virtualenv/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on master branch"

fi

# Report how big things are
echo "Docker cache space usage:"
du -hs $HOME/virtualenv/*
