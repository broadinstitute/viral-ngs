#!/bin/bash

set -e -o pipefail

mkdir -p /opt/viral-ngs/viral-ngs-etc
ln -s /opt/viral-ngs/source /opt/viral-ngs/viral-ngs-etc/viral-ngs
ln -s /opt/miniconda /opt/viral-ngs/viral-ngs-etc/conda-env
export VIRAL_CONDA_ENV_PATH=/opt/miniconda
mv /opt/docker/easy-deploy-viral-ngs.sh /opt/viral-ngs
mv /opt/docker/env_wrapper.sh /opt/viral-ngs

# setup/install viral-ngs
sync
/opt/viral-ngs/easy-deploy-viral-ngs.sh setup-git-local

# this not only prints the current version string, but it also saves it to the
# VERSION file for later use
PATH=/opt/miniconda/bin:$PATH
source activate /opt/miniconda
echo -n "viral-ngs version: "
/opt/viral-ngs/source/assembly.py --version
