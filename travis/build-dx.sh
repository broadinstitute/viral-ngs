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
      # blank this out until we're sure we want to test it this way...
      CMD_INPUT=""
    else
      CMD_INPUT=""
    fi

    defaults_json="pipes/WDL/dx-defaults-$workflow_name.json"
    if [ -f "$defaults_json" ]; then
      CMD_DEFAULTS="-defaults $defaults_json"
      # blank this out until we're sure this works...
      CMD_DEFAULTS=""
    else
      CMD_DEFAULTS=""
    fi

	  dx_id=$(java -jar dxWDL.jar compile \
      $workflow $CMD_INPUT $CMD_DEFAULTS -f \
      -destination /build/$VERSION/$workflow_name)
	  echo "Succeeded: $workflow_name = $dx_id"
    echo -e "$workflow_name\t$dx_id" >> $COMPILE_SUCCESS
  fi
done
# the presence of this file in the project denotes successful build
dx upload --brief --no-progress --destination /build/$VERSION/ $COMPILE_SUCCESS
