#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`

if [ "$GIT_BRANCH" != "master" ]; then
    echo "setting up Travis docker caches"
    mkdir -p caches/virtualenv caches/bin_bundles caches/tools_build

    rm -rf /home/travis/virtualenv
    ln caches/virtualenv /home/travis/virtualenv

    rm -rf bin_bundles
    ln caches/bin_bundles bin_bundles

    rm -rf tools/build
    ln caches/tools_build tools/build

else
    echo "Travis docker caches disabled for master branch"
fi
