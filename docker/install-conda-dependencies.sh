#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/viral-ngs),
# VIRAL_NGS_PATH (typically /opt/viral-ngs/source), and
# CONDA_DEFAULT_ENV (typically /opt/miniconda) to be set.
#
# A miniconda install must exist at $CONDA_DEFAULT_ENV
# and $CONDA_DEFAULT_ENV/bin must be in the PATH

set -e -o pipefail

echo "PATH:              ${PATH}"
echo "INSTALL_PATH:      ${INSTALL_PATH}"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "VIRAL_NGS_PATH:    ${VIRAL_NGS_PATH}"
echo "MINICONDA_PATH:    ${MINICONDA_PATH}"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV}"
CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"

# setup/install viral-ngs directory tree and conda dependencies
sync

REQUIREMENTS=""
for condafile in $*; do
	REQUIREMENTS="$REQUIREMENTS --file $condafile"
done

conda install -y -q $CONDA_CHANNEL_STRING -p "${CONDA_PREFIX}" $REQUIREMENTS

# clean up
conda clean -y --all
