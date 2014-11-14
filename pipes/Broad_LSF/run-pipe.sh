#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# load necessary Broad dotkits
eval `/broad/software/dotkit/init -b`
reuse -q LSF
reuse -q Python-3.4

# load Python virtual environment
source venv/bin/activate

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete \
	--jobs 100000 --immediate-submit \
	--jobscript bin/pipes/Broad_LSF/jobscript.sh \
	--cluster 'bin/pipes/Broad_LSF/cluster-submitter.py {dependencies} {config[logDir]}' \
	"$@"

