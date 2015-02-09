#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`

rm -rf tools/build
mkdir -p caches/tools_build

if [ "$GIT_BRANCH" != "master" ]; then
    echo "Travis docker caches allowed for branch $GIT_BRANCH"
    ln -s caches/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on master branch"
    mkdir -p tools/build

fi
