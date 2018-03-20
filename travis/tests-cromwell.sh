#!/bin/bash
set -x -e  # intentionally allow for pipe failures below

ln -s $GATK_PATH/GenomeAnalysisTK.jar .
mkdir -p workflows
cp *.jar test/input/WDL/test_inputs-*-local.json pipes/WDL/workflows/*.wdl pipes/WDL/workflows/tasks/*.wdl workflows
cd workflows

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test_inputs-$workflow_name-local.json"
	echo Looking at workflow $workflow_name
	echo dir is
	pwd
	if [ -f $input_json ]; then
		echo ${VIRAL_NGS_METADATA_PATH} > metadata_path.txt
		#export VIRAL_NGS_METADATA_VALUE_wdl_workflow=$workflow_name
		date
		echo "Executing $workflow_name using Cromwell on local instance"
		# the "cat" is to allow a pipe failure (otherwise it halts because of set -e)
		java -jar cromwell.jar run \
			$workflow_name.wdl \
			-i $input_json | tee cromwell.out
		if [ ${PIPESTATUS[0]} -gt 0 ]; then
			echo "error running $workflow_name"
			error_logs=$(grep stderr cromwell.out | perl -lape 's/.*\s(\S+)$/$1/g')
			for log in $error_logs; do
				echo "contents of stderr ($log):"
				cat $log | sed "s/^/[STDERR] /"
				echo "contents of stdout ($log):"
				cat `dirname $log`/stdout | sed "s/^/[STDOUT] /"
			done
			sync; sleep 30; exit 1
		fi
    fi
done

cd -
date
echo "note: there is no testing of output correctness yet..."
