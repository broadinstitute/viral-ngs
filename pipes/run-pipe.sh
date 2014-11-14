#!/bin/bash
# this is intended for use on the Broad LSF cluster

eval `/broad/software/dotkit/init -b`
reuse -q LSF
reuse -q Python-3.4

source venv/bin/activate

snakemake --timestamp --rerun-incomplete \
	--jobs 100000 --immediate-submit \
	--jobscript bin/pipes/jobscript.sh \
	--cluster 'bin/pipes/lsf-broad-submit.py {dependencies} {config[logDir]}' \
	"$@"

