#!/bin/sh
set -e

if [ ! -d $GATK_PATH -o ! -d $NOVOALIGN_PATH ]; then
    echo "Fetching encrypted Novoalign & GATK bundle for Travis"
    wget http://www.broadinstitute.org/~dpark/viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc
    openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in viral_ngs-gatk_novoalign-encrypted_for_travis.tar.gz.enc -out bin_bundles.tar.gz
    tar -xzpvf bin_bundles.tar.gz
fi

if [ $TRAVIS_PULL_REQUEST != "false" ]; then
    echo "This is a pull request: wiping out all installed tools and re-installing from scratch"
    rm -rf tools/build
    mkdir -p tools/build
fi

echo "Installing and validating bioinformatic tools"
python -m unittest -v test.test_tools.TestToolsInstallation
