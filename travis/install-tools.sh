#!/bin/bash
set -e

if [ ! -d $GATK_PATH ]; then
  if [ -z "$BUNDLE_SECRET" ]; then
    echo "ERROR: GATK is missing, but secret key is not set for auto-download."
    exit 1

  else
    echo "Fetching encrypted Novoalign & GATK bundle for Travis"
    pwd
    wget http://www.broadinstitute.org/~dpark/viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc
    openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc -out bin_bundles.tar.gz
    tar -xzpvf bin_bundles.tar.gz -C "$CACHE_DIR"

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
