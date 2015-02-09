#!/bin/bash
set -e

if [ ! -d $GATK_PATH -o ! -d $NOVOALIGN_PATH ]; then
    echo "Fetching encrypted Novoalign & GATK bundle for Travis"
    pwd
    wget http://www.broadinstitute.org/~dpark/viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc
    openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc -out bin_bundles.tar.gz
    tar -xzpvf bin_bundles.tar.gz
fi

echo "Installing and validating bioinformatic tools"
python -m unittest -v test.test_tools.TestToolsInstallation
