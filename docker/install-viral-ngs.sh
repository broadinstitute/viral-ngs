#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/viral-ngs)
# and VIRAL_NGS_PATH (typically /opt/viral-ngs/source) to be set.

set -e -o pipefail

export VIRAL_CONDA_ENV_PATH=/opt/miniconda

mkdir -p $INSTALL_PATH/viral-ngs-etc
ln -s $VIRAL_NGS_PATH $INSTALL_PATH/viral-ngs-etc/viral-ngs
ln -s $VIRAL_CONDA_ENV_PATH $INSTALL_PATH/viral-ngs-etc/conda-env
ln $VIRAL_NGS_PATH/easy-deploy-script/easy-deploy-viral-ngs.sh $INSTALL_PATH

# setup/install viral-ngs
sync
$INSTALL_PATH/easy-deploy-viral-ngs.sh setup-git-local

# this not only prints the current version string, but it also saves it to the
# VERSION file for later use
PATH=$VIRAL_CONDA_ENV_PATH/bin:$PATH
echo -n "viral-ngs version: "
$VIRAL_NGS_PATH/assembly.py --version
