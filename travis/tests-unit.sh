#!/bin/bash
#set -e

UNIT_TESTS_ALWAYS=1

if [ $UNIT_TESTS_ALWAYS ] || [ $TRAVIS_PULL_REQUEST == "false" ]; then
    export PYTEST_ADDOPTS="-rsxX -n 2 --durations=50 --fixture-durations=20 --junit-xml=pytest.xml --cov-report= --cov broad_utils --cov illumina --cov assembly --cov interhost --cov intrahost --cov metagenomics --cov ncbi --cov read_utils --cov reports --cov taxon_filter --cov tools --cov util --basetemp /dev/shm/test"
    pytest test/unit
else
    echo "Skipping unit tests for pull request - see branch tests."
fi
