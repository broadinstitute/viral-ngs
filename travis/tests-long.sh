#!/bin/bash
#set -e

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"
echo "TRAVIS_TAG: $TRAVIS_TAG"

if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then
    echo "This is on master or is a pull request: executing long running tests..."
    export PYTEST_ADDOPTS="-rsxX -n 2 --durations=50 --fixture-durations=20 --junit-xml=pytest.xml --cov-report= --cov broad_utils --cov illumina --cov assembly --cov interhost --cov intrahost --cov metagenomics --cov ncbi --cov read_utils --cov reports --cov taxon_filter --cov tools --cov util --basetemp /dev/shm/test"
    pytest --cov-append test/integration
    rc=$?; if [[ $rc != 0 ]]; then sleep 4 && exit $rc; fi
     # sleep to allow logs to be printed without truncation in the event of error
else
    echo "This is not a pull request: skipping long running tests."
fi
