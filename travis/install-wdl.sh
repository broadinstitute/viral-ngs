#!/bin/bash
set -e -o pipefail

echo "Fetching WDL related JAR files"
wget --quiet \
 https://github.com/dnanexus-rnd/dxWDL/releases/download/0.51/dxWDL-0.51.jar \
 https://github.com/broadinstitute/wdltool/releases/download/0.14/wdltool-0.14.jar \
 https://github.com/broadinstitute/cromwell/releases/download/29/cromwell-29.jar

echo "Fetching dx-toolkit"
TGZ=dx-toolkit-v0.240.1-ubuntu-14.04-amd64.tar.gz
wget --quiet https://wiki.dnanexus.com/images/files/$TGZ
tar -xzpf $TGZ
rm $TGZ

echo "Fetching quayctl"
wget --quiet wget https://github.com/coreos/quayctl/releases/download/v0.0.1/quayctl-linux-x64
mv quayctl-linux-64 quayctl
