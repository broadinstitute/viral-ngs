#!/bin/bash
set -e -o pipefail

GATK_VERSION=3.8
GATK_BUILD=1-0-gf15c1c3ef
GATK_URL="https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}.tar.bz2"
GATK_TAR=GenomeAnalysisTK-${GATK_VERSION}.tar
GATK_JAR=GenomeAnalysisTK.jar
GATK_JAR_MD5=186aee868bb7cffc18007966ace8d053

echo "Fetching GATK bundle"
pwd
wget --quiet "${GATK_URL}" -O ${GATK_TAR}

tar --warning=no-unknown-keyword -xvpf ${GATK_TAR}
rm ${GATK_TAR}
mv GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}/GenomeAnalysisTK.jar ${GATK_JAR}
rmdir GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}

md5sum --strict -c <(echo ${GATK_JAR_MD5} ${GATK_JAR})
gatk3-register ${GATK_JAR}
rm ${GATK_JAR}
