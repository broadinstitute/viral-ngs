#!/bin/bash
set -e

if [ $TRAVIS_PULL_REQUEST != "false" ]; then
    echo "This is a pull request: executing long running tests..."
    
    echo "(no long running tests currently exist... please add some)"
else
    echo "This is not a pull request: skipping long running tests."
fi
