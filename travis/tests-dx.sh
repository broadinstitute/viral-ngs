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

function dx_run_timeout_args {
    #
    # Construct command-line arguments for 'dx run' command
    # to set a timeout on the applets it runs
    #

    local dx_workflow_id="$1"
    local dx_extra_applet_id="$2"

    local dx_workflow_applet_ids=$(dx describe $dx_workflow_id | grep applet- | awk '{print $2;}')
    local dx_applet_ids="$dx_workflow_applet_ids $dx_extra_applet_id"
    local comma=""
    local timeout_args="{\"timeoutPolicyByExecutable\":{"
    for dx_applet_id in $dx_applet_ids
    do
        timeout_args="${timeout_args}${comma}\"$dx_applet_id\":{\"*\":{\"hours\":3}}"
        comma=","
    done
    timeout_args="$timeout_args}}"
    echo $timeout_args
}

TEST_LAUNCH_ALL="dxWDL-execute_all-launched.txt"
touch $TEST_LAUNCH_ALL
for workflow in pipes/WDL/workflows/*.wdl; do
  echo "testing $workflow..."
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
    input_json="test/input/WDL/test_inputs-$workflow_name-dnanexus.dx.json"
    if [ -f $input_json ]; then
       # launch simple test cases on DNAnexus CI project
       dx_workflow_id=$(grep -w "^$workflow_name" $COMPILE_SUCCESS | cut -f 2)
       timeout_args=$(dx_run_timeout_args $dx_workflow_id)
       dx_job_id=$(dx run \
           $dx_workflow_id -y --brief \
           -f $input_json \
           --name "$VERSION $workflow_name" \
           --destination /tests/$VERSION/$workflow_name \
           --extra-args $timeout_args \
           )
       if [ $? -eq 0 ]; then
           echo "Launched $workflow_name as $dx_job_id"
       else
           echo "Failed to build: $workflow_name"
       fi
       echo "Launched $workflow_name as $dx_job_id"
       echo -e "$workflow_name\t$dx_workflow_id\t$dx_job_id" >> $TEST_LAUNCH_ALL
    fi
  fi
done

# Special case: run test for the demux_launcher native applet (which invokes
# the demux_plus WDL workflow)
demux_launcher_id=$(grep demux_launcher $COMPILE_SUCCESS | cut -f 2)
demux_plus_workflow_id=$(grep demux_plus $COMPILE_SUCCESS | cut -f 2)
timeout_args=$(dx_run_timeout_args $demux_plus_workflow_id $demux_launcher_id)
dx_job_id=$(dx run \
  $demux_launcher_id -y --brief \
  -i upload_sentinel_record=record-Bv8qkgQ0jy198GK0QVz2PV8Y \
  --name "$VERSION demux_launcher" \
  -i folder=/tests/$VERSION/demux_launcher \
  --extra-args $timeout_args \
  )
echo "Launched demux_launcher as $dx_job_id"
echo -e "demux_launcher\t$demux_launcher_id\t$dx_job_id" >> $TEST_LAUNCH_ALL

# the presence of this file in the project denotes all tests launched
dx upload --brief --no-progress --destination /build/$VERSION/ $TEST_LAUNCH_ALL
