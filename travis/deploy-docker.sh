#!/bin/bash
set -e -o pipefail

if [ -n "$BUILD_PACKAGE" && -n "$DOCKER_PASS" && -n "$DOCKER_USER" ]; then
	echo "Deploying docker images to DockerHub"

    # deploy docker image to DockerHub
    echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER"
    if [ -n "$TRAVIS_TAG" ]; then
        DOCKER_REPO=broadinstitute/viral-ngs
        DOCKER_SHORT_TAG="latest"
        DOCKER_LONG_TAG="$TRAVIS_TAG"
    else
        if [ -z "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
            BRANCH_NAME="$TRAVIS_BRANCH"
        else
            BRANCH_NAME="$TRAVIS_PULL_REQUEST_BRANCH"
        fi
        DOCKER_REPO=broadinstitute/viral-ngs-dev
        DOCKER_SHORT_TAG="$(echo $TRAVIS_BRANCH | sed 's/-/_/g')"
        DOCKER_LONG_TAG="$(git describe --tags --always | sed 's/^v//' | perl -lape 's/(\d+.\d+.\d+)-/$1-dev-/')_$(echo $DOCKER_SHORT_TAG)"
    fi
    docker tag local/viral-ngs:build "$DOCKER_REPO:$DOCKER_SHORT_TAG"
    docker tag local/viral-ngs:build "$DOCKER_REPO:$DOCKER_LONG_TAG" 
    docker push "$DOCKER_REPO:$DOCKER_SHORT_TAG"
    docker push "$DOCKER_REPO:$DOCKER_LONG_TAG"

else
	echo "Skipping DockerHub deploy"

fi