#!/bin/sh
# properties = {properties}
# this is identical to the default jobscript with the exception of the exit code
{exec_job}

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
cat $1