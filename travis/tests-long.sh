#!/bin/sh

if [ $TRAVIS_PULL_REQUEST != "false" ]; then
    echo "This is a pull request, executing long running tests..."
    
    echo "(no long running tests currently exist... please add some)"
fi
