#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

# Always cache the encrypted bin_bundles that come from the Broad
# (pull requests from forks don't see the Travis secret key, so this
# allows external PRs to get tested properly)
mkdir -p $HOME/caches/bin_bundles
ln -s $HOME/caches/bin_bundles bin_bundles

# Only sometimes cache the tools/build directory
if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "Travis docker caches allowed for branch $TRAVIS_BRANCH"
    rm -rf tools/build
    mkdir -p $HOME/caches/tools_build
    ln -s $HOME/caches/tools_build tools/build

else
    echo "Travis docker cache disabled for tools/build on master branch"

fi

# Always report how big things are
echo "Docker cache space usage:"
du -hs $HOME/virtualenv $HOME/caches/*
