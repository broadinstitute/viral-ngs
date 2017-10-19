#!/bin/bash

echo "Script called to trigger tests in external repository..."

# parse out the second part of the job number (the '1' of '223.1')
JOB_NUMBER=$(echo $TRAVIS_JOB_NUMBER | perl -lape 's/^(\d+)\.(\d+)/$2/g')

# if this is the first job in the matrix, only then trigger the external 
# repository to avoid triple-building the Docker image
if [ "$JOB_NUMBER" == "1" ]; then
    # only initiate tests in other repo if the travis token string has a value
    if [ ! -z "$TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO" ]; then
        echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is defined. Triggering downstream repo..."
            # if this is a tagged release add that information to the dependent build request
        if [ -n "$TRAVIS_TAG" ]; then
            ./travis/trigger-travis.sh --script "env UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG ./travis/tests-unit.sh" broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO #"UPSTREAM_BRANCH=$TRAVIS_BRANCH UPSTREAM_TAG=$TRAVIS_TAG"
        else
            ./travis/trigger-travis.sh --script "env UPSTREAM_BRANCH=$TRAVIS_BRANCH ./travis/tests-unit.sh" broadinstitute viral-ngs-deploy $TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO #"UPSTREAM_BRANCH=$TRAVIS_BRANCH"
        fi
    else
        echo "TRAVIS_ACCESS_TOKEN_FOR_OTHER_REPO is undefined. Check the secure variable."
    fi
else
    echo "Not triggering, this is not the first job in the upstream build matrix."
fi
