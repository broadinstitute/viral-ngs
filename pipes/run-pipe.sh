#!/bin/bash

source venv/bin/activate

snakemake --timestamp --rerun-incomplete \
	--jobs 100000 --immediate-submit \
	--jobscript bin/pipes/jobscript.sh \
	--cluster 'bin/pipes/lsf-broad-submit.py {dependencies}' \
	"$@"

