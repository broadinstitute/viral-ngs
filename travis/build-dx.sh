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
dx select $DX_PROJECT

# compile with dxWDL
COMPILE_SUCCESS="dxWDL-compile_all-success.csv"
touch $COMPILE_SUCCESS
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
	  echo "Building $workflow to DNAnexus"
	  # TO DO: incorporate default file values once we figure out how
	  dx_id=`java -jar dxWDL-0.51.jar compile $workflow -destination /build/$VERSION/$workflow_name`
	  echo "Succeeded: $workflow_name = $dx_id"
    echo "$workflow_name,$dx_id" >> $COMPILE_SUCCESS
  fi
done
# the presence of this file in the project denotes successful build
dx upload --no-progress --destination /build/$VERSION/ $COMPILE_SUCCESS


## TO DO: trigger test executions on DNAnexus
#TEST_SUCCESS_ALL="dxWDL-execute_all-success.txt"
#for workflow in pipes/WDL/workflows/*.wdl; do
#  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
#    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
#  else
#    # launch simple test cases on DNAnexus CI project
#  fi
#done
## the presence of this file in the project denotes all tests pass
#dx upload --no-progress --destination /build/$VERSION/ $TEST_SUCCESS_ALL
