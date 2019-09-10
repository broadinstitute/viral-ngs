#!/bin/bash
set -e -o pipefail

if [[ -n "$DOCKER_PASS" && -n "$DOCKER_USER" && -n "$DOCKER_REGISTRY" ]]; then
	echo "Deploying docker image to $DOCKER_REGISTRY"

	# log in to Docker registry
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER" "$DOCKER_REGISTRY"

	if [ -n "$DOCKER_USER_MIRROR" -a -n "$DOCKER_PASS_MIRROR" -a -n "$DOCKER_REGISTRY_MIRROR" ]; then
		echo "$DOCKER_PASS_MIRROR" | docker login --password-stdin --username "$DOCKER_USER_MIRROR" "$DOCKER_REGISTRY_MIRROR"
	fi

	# tag and deploy
	for TAG in `travis/list-docker-tags.sh`; do
		echo "Pushing docker image to docker registry as $TAG"
		docker tag local/build-container:build $TAG
		docker push $TAG
	done

else
	echo "Skipping Docker image deploy unless DOCKER_REGISTRY, DOCKER_USER and DOCKER_PASS variables are set."
	exit 1

fi
