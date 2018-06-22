#!/bin/bash
set -e -o pipefail

starting_dir="$(pwd)"
test_dir="wdl_validate_test"

function cleanup(){
    echo "Cleaning up from WDL womtool; exit code: $?"
    cd "$starting_dir"
    if [ -d "$test_dir" ]; then
      rm -r "$test_dir"
    fi
}
trap cleanup EXIT SIGINT SIGQUIT SIGTERM

# ===== [ for WDL 1.0] =====
# validate each imported library of tasks on its own
# for tasks in pipes/WDL/workflows/tasks/*.wdl; do
#   echo "validating tasks $tasks"
#   if $(hash -r  womtool &> /dev/null); then
#     womtool validate $tasks
#   else
#     java -jar womtool.jar validate $tasks
#   fi
# done
# ==========================

mkdir "$test_dir"
cd "$test_dir"
cp ../pipes/WDL/workflows/tasks/*.wdl ../pipes/WDL/workflows/*.wdl .

# ===== [ for WDL draft-2 under womtool 32 ] =====
# only necessary for draft-2 WDL under womtool >30, with 1.0 task-only WDL files should validate
# the task validation block above can be uncommented and replace this with WDL 1.0
for workflow in ../pipes/WDL/workflows/tasks/*.wdl; do
  workflow=`basename $workflow`
  echo "validating $workflow"
  # include dummy workflow to validate under womtool 32
  printf "\n\nworkflow dummyworkflow_$(basename $workflow .wdl) {}" >> $workflow
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $workflow
  else
    java -jar ../womtool.jar validate $workflow
  fi
done
# ================================================

# copy in the original task.wdl files [remove with WDL 1.0] that lack dummy workflow
cp ../pipes/WDL/workflows/tasks/*.wdl ../pipes/WDL/workflows/*.wdl .
# validate the workflow files
# unfortunately, dxWDL now requires the -imports parameter and cromwell supports
# it as well but womtool validate does not yet support it! so we have to copy
# everything to a temp dir
for workflow in ../pipes/WDL/workflows/*.wdl; do
  workflow=`basename $workflow`
  echo "validating $workflow"
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $workflow
  else
    java -jar ../womtool.jar validate $workflow
  fi
done
