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
    tar -xzpvf GenomeAnalysisTK-3.6.tar.gz -C "$CACHE_DIR"

  fi
fi

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # encrypted bundle contains linux binary for Novoalign, remove that here
    unset NOVOALIGN_PATH
    # some conda packages dont exist on OSX
    cat requirements-conda.txt | grep -v diamond | grep -v kraken > $HOME/requirements-conda.txt
else
    # for linux, just use requirements-conda as-is
    cp requirements-conda.txt $HOME
fi

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default
conda create -y -m -c bioconda -p tools/conda-tools/default --file $HOME/requirements-conda.txt
./install_tools.py
