#!/bin/bash
set -e -o pipefail

# figure out appropriate tag names
if [ -n "$TRAVIS_TAG" ]; then
	# this is an official tagged release
	DOCKER_REPO="broadinstitute/viral-ngs"
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(echo $TRAVIS_TAG | sed 's/^v//')"
elif [ -n "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
	# this is a PR build (TRAVIS_BRANCH=master, TRAVIS_PULL_REQUEST_BRANCH=source of PR)
	DOCKER_REPO="broadinstitute/viral-ngs-dev"
	BRANCH_NAME="$TRAVIS_PULL_REQUEST_BRANCH"
	DOCKER_SHORT_TAG="$(echo $BRANCH_NAME | sed 's/-/_/g')-pull_request"
	DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g\S+/$1-beta$2/')_$(echo $BRANCH_NAME | sed 's/-/_/g')"
elif [[ "$TRAVIS_BRANCH" == "master" ]]; then
	# this is a master branch commit (e.g. merged pull request)
	DOCKER_REPO="broadinstitute/viral-ngs"
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g\S+/$1-rc$2/')"
else
	# this is an normal non-master branch commit
	DOCKER_REPO="broadinstitute/viral-ngs-dev"
	BRANCH_NAME="$TRAVIS_BRANCH"
	DOCKER_SHORT_TAG="$(echo $BRANCH_NAME | sed 's/-/_/g')"
	DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g\S+/$1-dev$2/')_$(echo $DOCKER_SHORT_TAG)"
fi

echo $DOCKER_REPO:$DOCKER_SHORT_TAG
echo $DOCKER_REPO:$DOCKER_LONG_TAG
