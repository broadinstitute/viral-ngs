#!/bin/bash
set -e -o pipefail

if [[ -n "$DOCKER_PASS" && -n "$DOCKER_USER" ]]; then
	echo "Deploying docker image to DockerHub"

	# log in to DockerHub
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER"

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
		#DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\S+)-(\d+)-g\S+/$1-beta-$2/')_$(echo $BRANCH_NAME | sed 's/-/_/g')"
	elif [[ "$TRAVIS_BRANCH" == "master" ]]; then
		# this is a master branch commit (e.g. merged pull request)
		DOCKER_REPO="broadinstitute/viral-ngs"
		DOCKER_SHORT_TAG="latest"
		DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\S+)-(\d+)-g\S+/$1-rc$2/')"
	else
		# this is an normal non-master branch commit
		DOCKER_REPO="broadinstitute/viral-ngs-dev"
		BRANCH_NAME="$TRAVIS_BRANCH"
		DOCKER_SHORT_TAG="$(echo $BRANCH_NAME | sed 's/-/_/g')"
		#DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1-dev-/')_$(echo $DOCKER_SHORT_TAG)"
	fi

	# tag and deploy
	for TAG in "$DOCKER_SHORT_TAG" "$DOCKER_LONG_TAG"; do
		if [ -n "$TAG" ]; then
			echo "Pushing docker image to DockerHub as $DOCKER_REPO:$TAG"
			docker tag local/viral-ngs:build "$DOCKER_REPO:$TAG" 
			docker push "$DOCKER_REPO:$TAG"
		fi
	done

else
	echo "Skipping DockerHub deploy unless DOCKER_USER and DOCKER_PASS variables are set."

fi
