#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# load necessary Broad dotkits
eval `/broad/software/dotkit/init -b`
reuse -q UGER
reuse -q Python-3.4
reuse -q Perl-5.10

# load config dirs from config.json
VENVDIR=`python -c 'import json;f=open("config.json");print(json.load(f)["venvDir"]);f.close()'`
BINDIR=`python -c 'import json;f=open("config.json");print(json.load(f)["binDir"]);f.close()'`

# load Python virtual environment
source "$VENVDIR/bin/activate"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
    --jobs 100000 --immediate-submit \
        --latency-wait 20 \
    --config mode=UGER \
    --directory . \
    --jobscript "$BINDIR/pipes/Broad_UGER/jobscript.sh" \
    --cluster $BINDIR'/pipes/Broad_UGER/cluster-submitter.py {dependencies} {config[logDir]}' \
    "$@"
