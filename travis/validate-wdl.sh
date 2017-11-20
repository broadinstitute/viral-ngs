#!/bin/bash
set -e -o pipefail

ln -s pipes/WDL/workflows/tasks .
for workflow in pipes/WDL/workflows/*.wdl; do
  java -jar wdltool-0.14.jar validate $workflow
done
rm tasks
