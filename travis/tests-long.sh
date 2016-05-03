#!/bin/bash
set -e

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"

if [ $TRAVIS_PULL_REQUEST != "false" -o $TRAVIS_BRANCH = "master" -o -n "$TRAVIS_TAG" ]; then
    echo "This is on master or is a pull request: executing long running tests..."
    py.test --cov-append test/integration
else
    echo "This is not a pull request: skipping long running tests."
fi
