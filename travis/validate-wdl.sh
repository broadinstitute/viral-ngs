#!/bin/bash
set -e -o pipefail

# validate each imported library of tasks on its own
for tasks in pipes/WDL/workflows/tasks/*.wdl; do
  echo "validating tasks $tasks"
  java -jar womtool.jar validate $tasks
done

# validate the workflow files
# unfortunately, dxWDL now requires the -imports parameter and cromwell supports
# it as well but womtool validate does not yet support it! so we have to copy
# everything to a temp dir
mkdir wdl_validate_test
cd wdl_validate_test
cp ../pipes/WDL/workflows/tasks/*.wdl ../pipes/WDL/workflows/*.wdl .
for workflow in ../pipes/WDL/workflows/*.wdl; do
  workflow=`basename $workflow`
  echo "validating $workflow"
  java -jar womtool.jar validate $workflow
done
cd -
rm -r wdl_validate_test
