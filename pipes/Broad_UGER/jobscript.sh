#!/bin/sh
# properties = {properties}
# this is identical to the default jobscript with the exception of the exit code

# Since UGER does not pass the environment cleanly, we have to set up dotkits here
eval `/broad/software/dotkit/init -b`
reuse -q UGER
reuse -q .python-3.4.3
reuse -q Perl-5.10
reuse -q Java-1.7
reuse -q .gcc-4.5.3
reuse -q .oracle-java-jdk-1.7.0-51-x86-64
reuse -q .bzip2-1.0.6 
reuse -q .zlib-1.2.6

{exec_job}

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null && exit 0 || exit 100