#!/bin/bash
set -e -o pipefail

NEW_TAG=`travis/list-docker-tags.sh | tail -1`
OLD_TAG="quay.io/broadinstitute/viral-ngs"

echo Replacing $OLD_TAG with $NEW_TAG in all task WDL files
sed -i -- "s|$OLD_TAG|$NEW_TAG|g" pipes/WDL/workflows/tasks/*.wdl

CURRENT_VERSION=`echo "$NEW_TAG" | cut -d\: -f2`
VERSION_PLACEHOLDER="viral-ngs_version_unknown"

echo Replacing $VERSION_PLACEHOLDER with $CURENT_VERSION in all task WDL files
sed -i -- "s|$VERSION_PLACEHOLDER|$CURRENT_VERSION|g" pipes/WDL/workflows/tasks/*.wdl

#echo "DEBUG: here are the modified runtime lines:"
#grep -n "$NEW_TAG" pipes/WDL/workflows/tasks/*.wdl
