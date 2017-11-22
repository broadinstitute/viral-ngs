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
COMPILE_SUCCESS="dxWDL-compile_all-success.txt"
touch $COMPILE_SUCCESS
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
	  echo "Building $workflow to DNAnexus"

    test_input_json_wdl="test/input/WDL/test_inputs-$workflow_name-dnanexus.json"
    if [ -f "$test_input_json_wdl" ]; then
      CMD_INPUT="-inputs $test_input_json_wdl"
      # blank this out until bugfix https://github.com/dnanexus-rnd/dxWDL/issues/69
      CMD_INPUT=""
    else
      CMD_INPUT=""
    fi

    defaults_json="pipes/WDL/dx-defaults-$workflow_name.json"
    if [ -f "$defaults_json" ]; then
      CMD_DEFAULTS="-defaults $defaults_json"
      # blank this out until bugfix https://github.com/dnanexus-rnd/dxWDL/issues/69
      CMD_DEFAULTS=""
    else
      CMD_DEFAULTS=""
    fi

	  dx_id=$(java -jar dxWDL-0.51.jar compile \
      $workflow $CMD_INPUT $CMD_DEFAULTS -f \
      -destination /build/$VERSION/$workflow_name)
	  echo "Succeeded: $workflow_name = $dx_id"
    echo -e "$workflow_name\t$dx_id" >> $COMPILE_SUCCESS
  fi
done
# the presence of this file in the project denotes successful build
dx upload --brief --no-progress --destination /build/$VERSION/ $COMPILE_SUCCESS


TEST_LAUNCH_ALL="dxWDL-execute_all-launched.txt"
touch $TEST_LAUNCH_ALL
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    input_json="test/input/WDL/test_inputs-$workflow_name-dnanexus.dx.json"
    if [ -f $input_json ]; then
       # launch simple test cases on DNAnexus CI project
       dx_workflow_id=$(grep $workflow_name $COMPILE_SUCCESS | cut -f 2)
       dx_job_id=$(dx run \
           $dx_workflow_id -y \
           -f $input_json \
           --name "$VERSION-$workflow_name" \
           --destination /tests/$VERSION/$workflow_name)
       echo "Launched $workflow_name as $dx_job_id"
       echo -e "$workflow_name\t$dx_workflow_id\t$dx_job_id" >> $TEST_LAUNCH_ALL
    fi
  fi
done
# the presence of this file in the project denotes all tests launched
dx upload --brief --no-progress --destination /build/$VERSION/ $TEST_LAUNCH_ALL
