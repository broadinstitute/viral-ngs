#!/bin/bash
set -e  # intentionally allow for pipe failures below

ln -s $GATK_PATH/GenomeAnalysisTK.jar .
mkdir -p workflows
ln -s pipes/WDL/workflows/*.wdl pipes/WDL/workflows/tasks/*.wdl workflows
cd workflows

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		input_json2="inputs-$workflow_name-local.json"
		sed s/METADATAPATH/${VIRAL_NGS_METADATA_PATH}/ > $input_json2
		echo PASSING INPUT:
		cat $input_json2
		#export VIRAL_NGS_METADATA_VALUE_wdl_workflow=$workflow_name
		date
		echo "Executing $workflow_name using Cromwell on local instance"
		# the "cat" is to allow a pipe failure (otherwise it halts because of set -e)
		java -jar cromwell.jar run \
			$workflow_name.wdl \
			-i $input_json2 | tee cromwell.out
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
