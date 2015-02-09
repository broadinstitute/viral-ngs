#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $HOME/caches/tools_build
    ln -s $HOME/caches/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on master branch"

fi
