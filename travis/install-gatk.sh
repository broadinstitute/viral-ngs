#!/bin/bash
set -e -o pipefail

GATK_VERSION=3.8
GATK_BUILD=1-0-gf15c1c3ef
GATK_URL="https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}.tar.bz2"
GATK_TAR=GenomeAnalysisTK-${GATK_VERSION}.tar
GATK_JAR=${CACHE_DIR}/GenomeAnalysisTK.jar
GATK_JAR_MD5=186aee868bb7cffc18007966ace8d053

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    MD5_PROG=md5
else
    MD5_PROG=md5sum
fi

if [ -f "${GATK_JAR}" ]; then
  echo "GATK already exists in cache, skipping download"

else
  echo "Fetching GATK bundle for Travis"
  pwd
  wget --quiet "${GATK_URL}" -O ${GATK_TAR}

  # It appears that GATK tarball is produced on OS X leading to warnings
  TAR_OPTS=
  [[ "$TRAVIS_OS_NAME" = "linux" ]] && TAR_OPTS="--warning=no-unknown-keyword"
  tar "$TAR_OPTS" -xvpf ${GATK_TAR} -C "$CACHE_DIR"
  mv ${CACHE_DIR}/GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}/GenomeAnalysisTK.jar ${GATK_JAR}
  rmdir ${CACHE_DIR}/GenomeAnalysisTK-${GATK_VERSION}-${GATK_BUILD}
fi

${MD5_PROG} --strict -c <(echo ${GATK_JAR_MD5} ${GATK_JAR})
