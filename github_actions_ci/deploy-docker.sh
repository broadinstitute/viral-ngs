#!/bin/bash
set -e -x -o pipefail

echo "Deploying docker image to $DOCKER_REGISTRY"
# tag and deploy
for TAG in `github_actions_ci/list-docker-tags.sh`; do
	echo "Pushing docker image to docker registry as $TAG"
	docker tag local/build-container:build $TAG
	docker push $TAG # ToDo: uncomment when tag is parsed correctly.
done