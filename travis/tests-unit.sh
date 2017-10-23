#!/bin/bash

UNIT_TESTS_ALWAYS=1

if [ $PERFORM_TESTING = "true" ]; then
    if [ $UNIT_TESTS_ALWAYS ] || [ $TRAVIS_PULL_REQUEST = "false" ]; then
        export PYTEST_ADDOPTS="-rsxX -n 2 --durations=50 --fixture-durations=20 --junit-xml=pytest.xml --cov-report= --cov broad_utils --cov illumina --cov assembly --cov interhost --cov intrahost --cov metagenomics --cov ncbi --cov read_utils --cov reports --cov taxon_filter --cov tools --cov util --basetemp /dev/shm/test"

        pytest test/unit
        rc=$?; if [[ $rc != 0 ]]; then sleep 4 && exit $rc; fi
        # sleep to allow logs to be printed without truncation in the event of error
    else
        echo "Skipping unit tests for pull request - see branch tests."
    fi
else
    echo "Unit tests not performed for this slot of the build matrix."
fi