#!/bin/bash

echo "Script called to trigger tests in external repository..."

# only initiate tests in other repo if the travis token string has a value
if [ ! -z "$TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO" ]; then
    echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is defined"
    # if this is on mater or a pull request, but not a tagged release
    if [ $TRAVIS_PULL_REQUEST != "false" -o $TRAVIS_BRANCH = "master" -o -n "$TRAVIS_TAG" ]; then
        echo "This is on master or is a pull request: executing long running tests..."
        # if this is a tagged release add that information to the dependent build request
        if [ -n "$TRAVIS_TAG" ]; then
            ./travis/trigger-travis.sh --script 'UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG sh travis/tests-long.sh' broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO "UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG"
        else
            ./travis/trigger-travis.sh --script 'UPSTREAM_BRANCH=$TRAVIS_BRANCH sh travis/tests-long.sh' broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO "UPSTREAM_BRANCH=$TRAVIS_BRANCH"
        fi
    fi
else
    echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is undefined. Check the secure variable."
fi