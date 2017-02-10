#!/bin/bash
set -e

if [ ! -d $GATK_PATH ]; then
  if [ -z "$BUNDLE_SECRET" ]; then
    echo "ERROR: GATK is missing, but secret key is not set for auto-download."
    exit 1

  else
    echo "Fetching encrypted GATK bundle for Travis"
    pwd
    wget https://storage.googleapis.com/sabeti-public/software_testing/GenomeAnalysisTK-3.6.tar.gz.enc
    openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in GenomeAnalysisTK-3.6.tar.gz.enc -out GenomeAnalysisTK-3.6.tar.gz
    md5sum GenomeAnalysisTK-3.6.tar.gz
    # It appears that GATK tarball is produced on OS X leading to warnings
    TAR_OPTS=
    [[ "$TRAVIS_OS_NAME" = "linux" ]] && TAR_OPTS="--warning=no-unknown-keyword"
    tar "$TAR_OPTS" -xzpvf GenomeAnalysisTK-3.6.tar.gz -C "$CACHE_DIR"

  fi
fi


if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # encrypted bundle contains linux binary for Novoalign, remove that here
    unset NOVOALIGN_PATH
    # some conda packages dont exist on OSX
    cat requirements-conda.txt | grep -v kraken > $HOME/requirements-conda.txt
else
    # for linux, just use requirements-conda as-is
    cp requirements-conda.txt $HOME
fi

# Set to conda's java
export JAVA_HOME="$(pwd)/tools/conda-tools/default/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

for i in $(seq 3); do
  conda create -y -m -c bioconda -p tools/conda-tools/default --file $HOME/requirements-conda.txt python="$TRAVIS_PYTHON_VERSION" && break
  sleep 5
done

echo 'Sourcing default environment'
source activate tools/conda-tools/default
conda info -a # for debugging

./install_tools.py
