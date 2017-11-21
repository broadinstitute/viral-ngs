#!/bin/bash
set -e -o pipefail

if [[ -n "$DOCKER_PASS" && -n "$DOCKER_USER" ]]; then
	echo "Deploying docker image to DockerHub"

	# log in to DockerHub
	echo "$DOCKER_PASS" | docker login --password-stdin --username "$DOCKER_USER"

	# tag and deploy
	FIRST_TIME=1
	for TAG in `travis/list-docker-tags.sh`; do
		echo "Pushing docker image to DockerHub as $TAG"
		docker tag local/viral-ngs:build $TAG
		if [ $FIRST_TIME -gt 0 ]; then
			FIRST_TIME=0
			docker push $TAG
		else
			# this should be instantaneous anyway, and doesn't penalize our travis log limit
			docker push $TAG > /dev/null
		fi
	done

else
	echo "Skipping DockerHub deploy unless DOCKER_USER and DOCKER_PASS variables are set."

fi
