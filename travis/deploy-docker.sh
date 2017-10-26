#!/bin/bash
set -e -o pipefail

if [[ -n "$BUILD_PACKAGE" && -n "$DOCKER_PASS" && -n "$DOCKER_USER" ]]; then
	echo "Deploying docker image to DockerHub"

	# log in to DockerHub
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER"

	# figure out appropriate tag names
	if [ -n "$TRAVIS_TAG" ]; then
		# this is an official tagged release
		DOCKER_REPO="broadinstitute/viral-ngs"
		DOCKER_SHORT_TAG="latest"
		DOCKER_LONG_TAG="$TRAVIS_TAG"
	else
		DOCKER_REPO="broadinstitute/viral-ngs-dev"
		if [ -n "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
			# this is a PR build
			BRANCH_NAME="$TRAVIS_PULL_REQUEST_BRANCH"
			DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1-beta-/')_$(echo $BRANCH_NAME | sed 's/-/_/g')"
		else
			# this is an normal branch commit
			BRANCH_NAME="$TRAVIS_BRANCH"
			DOCKER_SHORT_TAG="$(echo $BRANCH_NAME | sed 's/-/_/g')"
			if [[ "$BRANCH_NAME" == "master" ]]; then
				DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1-rc-/')_$(echo $DOCKER_SHORT_TAG)"
			else
				DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1-dev-/')_$(echo $DOCKER_SHORT_TAG)"
			fi
		fi
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
	echo "Skipping DockerHub deploy"

fi
