#!/bin/bash
set -e -o pipefail

if [ -z "$BUNDLE_SECRET" ]; then
  echo "ERROR: GATK is missing, but secret key is not set for auto-download."
  exit 1

elif [ -f "$CACHE_DIR/GenomeAnalysisTK.jar" ]; then
  echo "GATK already exists in cache, skipping download"

else
  echo "Fetching encrypted GATK bundle for Travis"
  pwd
  wget --quiet https://storage.googleapis.com/sabeti-public/software_testing/GenomeAnalysisTK-3.6.tar.gz.enc
  openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in GenomeAnalysisTK-3.6.tar.gz.enc -out GenomeAnalysisTK-3.6.tar.gz
  if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    md5 GenomeAnalysisTK-3.6.tar.gz
  else
    md5sum GenomeAnalysisTK-3.6.tar.gz
  fi
  # It appears that GATK tarball is produced on OS X leading to warnings
  TAR_OPTS=
  [[ "$TRAVIS_OS_NAME" = "linux" ]] && TAR_OPTS="--warning=no-unknown-keyword"
  tar "$TAR_OPTS" -xzvpf GenomeAnalysisTK-3.6.tar.gz -C "$CACHE_DIR"

fi
