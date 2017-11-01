#!/bin/bash

set -e -o pipefail

if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then
    pushd docs
    make html && echo "Docs built successfully!" || echo "Docs did NOT build successfully."
    popd
fi