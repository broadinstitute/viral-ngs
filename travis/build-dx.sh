#!/bin/bash
set -e -o pipefail

if [ -z "$DX_API_TOKEN" ]; then
  echo "ERROR: DX_API_TOKEN is not set, this is needed to build dxWDL workflows."
  exit 1
fi

# obtain version tag
VERSION=`travis/list-docker-tags.sh | tail -1 | sed 's/:/\//'`

# log in to DNAnexus
source dx-toolkit/environment
dx login --token "$DX_API_TOKEN" --noprojects
dx select project-F856jv809y3VkzyFGkKqX367

# compile with dxWDL
for workflow in pipes/WDL/workflows/*.wdl; do
  # TO DO: incorporate default file values once we figure out how
  java -jar dxWDL-0.51.jar compile $workflow -destination /build/$VERSION
done

# TO DO: trigger test executions on DNAnexus
