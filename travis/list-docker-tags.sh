#!/bin/bash
set -e -o pipefail

if [ -z "$DOCKER_REPO_PROD" -o -z "$DOCKER_REPO_DEV" -o -z "$TRAVIS_BRANCH" ]; then
	echo "This script requires the following environment variables set: DOCKER_REPO_PROD, DOCKER_REPO_DEV, and TRAVIS_BRANCH."
	exit 1
fi

# figure out appropriate tag names
if [ -n "$TRAVIS_TAG" ]; then
	# this is an official tagged release
	DOCKER_REPO=$DOCKER_REPO_PROD
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(echo $TRAVIS_TAG | sed 's/^v//')"
elif [ -n "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
	# this is a PR build (TRAVIS_BRANCH=target of PR, TRAVIS_PULL_REQUEST_BRANCH=source of PR)
	DOCKER_REPO=$DOCKER_REPO_DEV
	BRANCH_NAME="$TRAVIS_PULL_REQUEST_BRANCH"
	DOCKER_SHORT_TAG="$BRANCH_NAME-pull_request"
	DOCKER_LONG_TAG="$(git describe --tags --always | sed s/^v//)-$(echo $DOCKER_SHORT_TAG)"
	#DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g(\S+)/$1-beta$2-g$3/')-$(echo $BRANCH_NAME | sed 's/-/_/g')"
elif [[ "$TRAVIS_BRANCH" == "master" ]]; then
	# this is a master branch commit (e.g. merged pull request)
	DOCKER_REPO=$DOCKER_REPO_PROD
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g\S+/$1-rc$2/')"
else
	# this is an normal non-master branch commit
	DOCKER_REPO=$DOCKER_REPO_DEV
	BRANCH_NAME="$TRAVIS_BRANCH"
	DOCKER_SHORT_TAG="$BRANCH_NAME"
	DOCKER_LONG_TAG="$(git describe --tags --always | sed s/^v//)-$(echo $DOCKER_SHORT_TAG)"
	#DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g(\S+)/$1-dev$2-g$3/')-$(echo $DOCKER_SHORT_TAG)"
fi

echo $DOCKER_REPO:$DOCKER_SHORT_TAG
echo $DOCKER_REPO:$DOCKER_LONG_TAG
