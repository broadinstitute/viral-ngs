#!/bin/bash
set -e -o pipefail

NEW_TAG=`travis/list-docker-tags.sh | tail -1`
OLD_TAG="broadinstitute/viral-ngs"

echo Replacing $OLD_TAG with $NEW_TAG in all task WDL files
sed -i -- "s|$OLD_TAG|$NEW_TAG|g" pipes/WDL/workflows/tasks/*.wdl

echo "DEBUG: here are the modified runtime lines:"
grep \'$NEW_TAG\' pipes/WDL/workflows/tasks/*.wdl
