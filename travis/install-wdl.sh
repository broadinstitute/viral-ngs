#!/bin/bash
set -e -o pipefail

if [ -z "$BUNDLE_SECRET" ]; then
  echo "ERROR: GATK is missing, but secret key is not set for auto-download."
  exit 1
fi

echo "Fetching WDL related JAR files"
wget --quiet \
 https://github.com/dnanexus-rnd/dxWDL/releases/download/0.51/dxWDL-0.51.jar \
 https://github.com/broadinstitute/wdltool/releases/download/0.14/wdltool-0.14.jar \
 https://github.com/broadinstitute/cromwell/releases/download/29/cromwell-29.jar \
 https://storage.googleapis.com/sabeti-public/software_testing/GenomeAnalysisTK-3.6.tar.gz.enc

echo "Unpacking GATK"
openssl aes-256-cbc -d -k "$BUNDLE_SECRET" -in GenomeAnalysisTK-3.6.tar.gz.enc -out GenomeAnalysisTK-3.6.tar.gz
# It appears that GATK tarball is produced on OS X leading to warnings
tar --warning=no-unknown-keyword -xzpf GenomeAnalysisTK-3.6.tar.gz

echo "Fetching dx-toolkit"
TGZ=dx-toolkit-v0.240.1-ubuntu-14.04-amd64.tar.gz
wget --quiet https://wiki.dnanexus.com/images/files/$TGZ
tar -xzpf $TGZ
rm $TGZ
