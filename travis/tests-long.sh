#!/bin/bash
set -e

echo travis_fold:start:tests-long

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"

if [ $TRAVIS_PULL_REQUEST != "false" -o $TRAVIS_BRANCH = "master" -o -n "$TRAVIS_TAG" ]; then
    echo "This is on master or is a pull request: executing long running tests..."
    nosetests -v --with-xunit --with-coverage --nocapture \
        --cover-inclusive --cover-branches --cover-tests \
        --cover-package broad_utils,illumina,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util \
        -w test/integration/

else
    echo "This is not a pull request: skipping long running tests."
fi

echo travis_fold:end:tests-long
