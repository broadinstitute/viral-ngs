#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# load necessary Broad dotkits
eval `/broad/software/dotkit/init -b`
reuse -q UGER
#reuse -q Python-3.4
reuse -q .python-3.4.3
reuse -q Perl-5.10
reuse -q Java-1.7
reuse -q .gcc-4.5.3
reuse -q .oracle-java-jdk-1.7.0-51-x86-64
reuse -q .bzip2-1.0.6 
reuse -q .zlib-1.2.6


# load config dirs from config.json
VENVDIR=`python -c 'import json; import os; f=open("config.json");print(os.path.realpath(json.load(f)["venvDir"]));f.close()'`
BINDIR=`python -c 'import json; import os; f=open("config.json");print(os.path.realpath(json.load(f)["binDir"]));f.close()'`

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
