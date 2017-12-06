#!/bin/bash
set -e -o pipefail

# obtain version tag
VERSION=`travis/list-docker-tags.sh | tail -1 | sed 's/:/\//'`

# log in to DNAnexus
source dx-toolkit/environment
if [ -n "$DX_API_TOKEN" ]; then
  dx login --token "$DX_API_TOKEN" --noprojects
  dx select $DX_PROJECT
fi

COMPILE_SUCCESS="dxWDL-compile_all-success.txt"
if [ ! -f $COMPILE_SUCCESS ]; then
  dx download --no-progress /build/$VERSION/$COMPILE_SUCCESS
fi

TEST_LAUNCH_ALL="dxWDL-execute_all-launched.txt"
touch $TEST_LAUNCH_ALL
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
    input_json="test/input/WDL/test_inputs-$workflow_name-dnanexus.dx.json"
    if [ -f $input_json ]; then
       # launch simple test cases on DNAnexus CI project
       dx_workflow_id=$(grep $workflow_name $COMPILE_SUCCESS | cut -f 2)
       dx_job_id=$(dx run \
           $dx_workflow_id -y --brief \
           -f $input_json \
           --name "$VERSION $workflow_name" \
           --destination /tests/$VERSION/$workflow_name)
       echo "Launched $workflow_name as $dx_job_id"
       echo -e "$workflow_name\t$dx_workflow_id\t$dx_job_id" >> $TEST_LAUNCH_ALL
    fi
  fi
done
# the presence of this file in the project denotes all tests launched
dx upload --brief --no-progress --destination /build/$VERSION/ $TEST_LAUNCH_ALL
