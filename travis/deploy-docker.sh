#!/bin/bash
set -e -o pipefail

if [[ -n "$DOCKER_PASS" && -n "$DOCKER_USER" ]]; then
	echo "Deploying docker image to DockerHub"

	# log in to DockerHub
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER"

	# tag and deploy
	for TAG in `travis/list-docker-tags.sh`; do
		echo "Pushing docker image to DockerHub as $TAG"
		docker tag local/viral-ngs:build $TAG 
		docker push $TAG
	done

else
	echo "Skipping DockerHub deploy unless DOCKER_USER and DOCKER_PASS variables are set."

fi
