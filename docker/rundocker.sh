#!/usr/bin/env bash

# A wrapper script to run viral-ngs docker images
# The following paths have to be modified according to end-user environment
NOVOALIGN_PATH="/opt/novocraft" # Directory where novoalign.lic license file can befound
GATK_PATH="/opt/GenomeAnalysisTK-3.8" # Directory where the correct GATK jar file can be found
IMAGE_HASH_OR_TAG="quay.io/broadinstitute/viral-ngs:latest" # This can be found by running this command 'docker images'
DATA_DIR="$1"; shift
GID=$(id -g $USER)

docker run -it \
       --volume $NOVOALIGN_PATH:/novoalign \
       --volume $GATK_PATH:/gatk \
       --volume $DATA_DIR:/user-data \
       -e RUN_USER_ID=$UID -e RUN_GROUP_ID=$GID \
       $IMAGE_HASH_OR_TAG \
       "$@"
