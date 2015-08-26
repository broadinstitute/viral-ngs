#!/bin/sh
# properties = {properties}
# this is identical to the default jobscript with the exception of the exit code

BINDIR=`python -c 'import json; import os; f=open("config.json");print(os.path.realpath(json.load(f)["binDir"]));f.close()'`

source "$BINDIR/pipes/Broad_UGER/setup_dotkits.sh"

{exec_job}

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null && exit 0 || exit 100