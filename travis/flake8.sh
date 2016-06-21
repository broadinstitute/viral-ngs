#!/bin/bash
set -e

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"

if [ $TRAVIS_PULL_REQUEST != "false" -o $TRAVIS_BRANCH = "master" -o -n "$TRAVIS_TAG" ]; then
    echo "This is on master or is a pull request: checking code style..."
    pip install 'flake8<=3'
    flake8 --exit-zero .
else
    echo "This is not a pull request: skipping code style check."
fi
