#!/bin/bash
set -e -o pipefail

ln -s pipes/WDL/workflows pipes/WDL/workflows/tasks .

for workflow in pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_values-$workflow_name-local.json"
	if [ -f $input_json ]; then
		java -jar cromwell-29.jar run workflows/$workflow_name.wdl -i $input_json

		echo "note: there is no testing of output yet..."
    fi
done
