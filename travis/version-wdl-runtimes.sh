#!/bin/bash
set -e -o pipefail

NEW_TAG=`travis/list-docker-tags.sh | tail -1`
OLD_STRING='docker: "broadinstitute/viral-ngs"'
NEW_STRING='docker: "'$NEW_TAG'"'

echo Replacing $OLD_STRING with $NEW_STRING in all task WDL files
sed -i -- "s/$OLD_STRING/$NEW_STRING/g" pipes/workflows/tasks/*.wdl

echo "DEBUG: here are the modified runtime lines:"
grep "docker:" pipes/workflows/tasks/*.wdl
