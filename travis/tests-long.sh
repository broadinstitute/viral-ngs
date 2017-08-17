#!/bin/bash
#set -e

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"
echo "TRAVIS_TAG: $TRAVIS_TAG"

if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then
    echo "This is on master or is a pull request: executing long running tests..."
    # pytest uses "-r[outputLetterCodes]" to set outputs, these are not other arguments composed into a single argument
    export PYTEST_ADDOPTS="-rsxX --exitfirst --durations=50 --fixture-durations=20 --junit-xml=pytest.xml --cov-report= --cov broad_utils --cov illumina --cov assembly --cov interhost --cov intrahost --cov metagenomics --cov ncbi --cov read_utils --cov reports --cov taxon_filter --cov tools --cov util"
    pytest --cov-append test/integration || exit $?
    #sleep 4 # sleep to allow logs to be printed without truncation in the event of error
else
    echo "This is not a pull request: skipping long running tests."
fi
