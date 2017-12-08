#!/bin/bash
set -e -o pipefail

if [[ -n "$DOCKER_PASS" && -n "$DOCKER_USER" && -n "$DOCKER_REGISTRY" ]]; then
	echo "Deploying docker image to $DOCKER_REGISTRY"

	# log in to Docker registry
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER" "$DOCKER_REGISTRY"
	echo "$DOCKER_PASS_MIRROR" | docker login --password-stdin --username "$DOCKER_USER_MIRROR" "$DOCKER_REGISTRY_MIRROR"

	# tag and deploy
	for TAG in `travis/list-docker-tags.sh`; do
		echo "Pushing docker image to docker registry as $TAG"
		docker tag local/viral-ngs:build $TAG
		docker push $TAG
	done

else
	echo "Skipping Docker image deploy unless DOCKER_REGISTRY, DOCKER_USER and DOCKER_PASS variables are set."

fi
