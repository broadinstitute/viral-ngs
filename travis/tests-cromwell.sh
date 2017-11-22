#!/bin/bash
set -e  # intentionally allow for pipe failures below

ln -s pipes/WDL/workflows pipes/WDL/workflows/tasks .

for workflow in pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "Executing $workflow_name using Cromwell on local instance"
		# the "cat" is to allow a pipe failure (otherwise it halts because of set -e)
		java -jar cromwell-29.jar run \
			workflows/$workflow_name.wdl \
			-i $input_json | cat > cromwell.out
		error_logs=$(grep stderr cromwell.out | perl -lape 's/.*\s(\S+)$/$1/g')
		if [ -n "$error_logs" ]; then
			echo "error running $workflow_name"
			for log in $error_logs; do
				echo "contents of stderr ($log):" >&2
				cat $log | sed "s/^/[STDERR] /" >&2
			done
			sync; sleep 30; exit 1
		fi
    fi
done

date
echo "note: there is no testing of output correctness yet..."
