#!/bin/bash
set -e -o pipefail

starting_dir="$(pwd)"
test_dir="wdl_validate_test"

function cleanup(){
    echo "Cleaning up from WDL womtool; exit code: $?"
    cd "$starting_dir"
    rm -r $test_dir
}
trap cleanup EXIT SIGINT SIGQUIT SIGTERM

# validate each imported library of tasks on its own
for tasks in pipes/WDL/workflows/tasks/*.wdl; do
  echo "validating tasks $tasks"
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $tasks
  else
    java -jar womtool.jar validate $tasks
  fi
done

# validate the workflow files
# unfortunately, dxWDL now requires the -imports parameter and cromwell supports
# it as well but womtool validate does not yet support it! so we have to copy
# everything to a temp dir
mkdir "$test_dir"
cd "$test_dir"
cp ../pipes/WDL/workflows/tasks/*.wdl ../pipes/WDL/workflows/*.wdl .
for workflow in ../pipes/WDL/workflows/*.wdl; do
  workflow=`basename $workflow`
  echo "validating $workflow"
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $workflow
  else
    java -jar ../womtool.jar validate $workflow
  fi
done
