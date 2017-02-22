#!/bin/bash

echo "Script called to trigger tests in external repository..."

# only initiate tests in other repo if the travis token string has a value
if [ ! -z "$TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO" ]; then
    echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is defined"
        # if this is a tagged release add that information to the dependent build request
    if [ -n "$TRAVIS_TAG" ]; then
        ./travis/trigger-travis.sh --script "env UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG bash ./travis/tests-unit.sh" broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO #"UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG"
    else
        ./travis/trigger-travis.sh --script "env UPSTREAM_BRANCH=$TRAVIS_BRANCH bash ./travis/tests-unit.sh" broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO #"UPSTREAM_BRANCH=$TRAVIS_BRANCH"
    fi
else
    echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is undefined. Check the secure variable."
fi