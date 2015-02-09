#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

rm -rf tools/build
mkdir -p caches/tools_build

if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    ln -s caches/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on master branch"
    mkdir -p tools/build

fi
