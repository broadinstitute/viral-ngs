#!/bin/bash
set -e -o pipefail

for workflow in pipes/WDL/workflows/*.wdl; do
  java -jar wdltool-0.14.jar validate $workflow
done
