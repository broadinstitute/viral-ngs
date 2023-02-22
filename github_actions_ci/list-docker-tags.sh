#!/bin/bash
#
# This script lists all docker tags that should apply to this particular
# deploy, including short tags (like repo:latest) and long tags (like
# repo:full-version-number). Tags are emitted one per line.
#
# Important note: the last tag is considered the definitive, uniquely
# identifiable tag that should be used for all build steps that depend on
# the docker image (e.g. version-wdl-runtimes.sh). All others are aliases.

set -e -o pipefail

if [ -z "$DOCKER_REPO_PROD" -o -z "$DOCKER_REPO_DEV" -o -z "$GITHUB_ACTIONS_BRANCH" ]; then
	echo "This script requires the following environment variables set: DOCKER_REPO_PROD, DOCKER_REPO_DEV, and GITHUB_ACTIONS_BRANCH."
	exit 1
fi

# figure out appropriate tag names
if [ -n "$GITHUB_ACTIONS_TAG" ]; then
	# this is an official tagged release
	DOCKER_REPO=$DOCKER_REPO_PROD
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(echo $GITHUB_ACTIONS_TAG | sed 's/^v//')"
	PUSH_TO_MIRROR=true
elif [ -n "$GITHUB_ACTIONS_PULL_REQUEST_BRANCH" ]; then
	# this is a PR build (GITHUB_ACTIONS_BRANCH=target of PR, GITHUB_ACTIONS_PULL_REQUEST_BRANCH=source of PR)
	DOCKER_REPO=$DOCKER_REPO_DEV
	BRANCH_NAME="$GITHUB_ACTIONS_PULL_REQUEST_BRANCH"
	DOCKER_SHORT_TAG="$BRANCH_NAME-pull_request"
	DOCKER_LONG_TAG="$(git describe --tags --always | sed s/^v//)-$(echo $DOCKER_SHORT_TAG)"
	#DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g(\S+)/$1-beta$2-g$3/')-$(echo $BRANCH_NAME | sed 's/-/_/g')"
elif [[ "$GITHUB_ACTIONS_BRANCH" == "master" ]]; then
	# this is a master branch commit (e.g. merged pull request)
	DOCKER_REPO=$DOCKER_REPO_PROD
	DOCKER_SHORT_TAG="latest"
	DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g\S+/$1-rc$2/')"
	PUSH_TO_MIRROR=true
else
	# this is an normal non-master branch commit
	DOCKER_REPO=$DOCKER_REPO_DEV
	BRANCH_NAME="$GITHUB_ACTIONS_BRANCH"
	DOCKER_SHORT_TAG="$BRANCH_NAME"
	DOCKER_LONG_TAG="$(git describe --tags --always | sed s/^v//)-$(echo $DOCKER_SHORT_TAG)"
	#DOCKER_LONG_TAG="$(git describe --tags --always | perl -lape 's/^v?(\S+)-(\d+)-g(\S+)/$1-dev$2-g$3/')-$(echo $DOCKER_SHORT_TAG)"
fi

if [ -n "$DOCKER_REPO_PROD_MIRROR" -a -n "$PUSH_TO_MIRROR" ]; then
	# if this should be sent to prod, also echo the mirror registry tags
	echo $DOCKER_REPO_PROD_MIRROR:$DOCKER_SHORT_TAG
	echo $DOCKER_REPO_PROD_MIRROR:$DOCKER_LONG_TAG
fi

echo $DOCKER_REPO:$DOCKER_SHORT_TAG
echo $DOCKER_REPO:$DOCKER_LONG_TAG
