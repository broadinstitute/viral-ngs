#!/bin/bash
# Wrappers around Snakemake for use on the Broad LSF cluster

# load necessary Broad dotkits
eval `/broad/software/dotkit/init -b`
reuse -q LSF
reuse -q Python-3.4
reuse -q Perl-5.10

# load config dirs from config.yaml
VENVDIR=`python -c 'import yaml;f=open("config.yaml");print(yaml.safe_load(f)["venv_dir"]);f.close()'`
BINDIR=`python -c 'import yaml;f=open("config.yaml");print(yaml.safe_load(f)["bin_dir"]);f.close()'`
DATADIR=`python -c 'import yaml; import os; f=open("config.yaml");print(os.path.realpath(yaml.safe_load(f)["data_dir"]));f.close()'`

# load Python virtual environment
source "$VENVDIR/bin/activate"

# invoke Snakemake in cluster mode with custom wrapper scripts
snakemake --timestamp --rerun-incomplete --keep-going --nolock \
	--jobs 100000 --immediate-submit \
        --latency-wait 20 \
	--config mode=LSF job_profiler="$BINDIR/pipes/Broad_LSF/lsf-report.py" \
	--directory . \
	--jobscript "$BINDIR/pipes/Broad_LSF/jobscript.sh" \
	--cluster $BINDIR"/pipes/Broad_LSF/cluster-submitter.py {dependencies} $DATADIR {config[log_dir]}" \
	"$@"
