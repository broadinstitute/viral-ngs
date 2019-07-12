#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/viral-ngs),
# VIRAL_NGS_PATH (typically /opt/viral-ngs/source), and
# CONDA_DEFAULT_ENV (typically /opt/miniconda) to be set.
#
# A miniconda install must exist at $CONDA_DEFAULT_ENV
# and $CONDA_DEFAULT_ENV/bin must be in the PATH
#
# Otherwise, this only requires the existence of the following files:
#	requirements-minimal.txt
#	requirements-conda.txt
#	requirements-conda-tests.txt
#	requirements-py3.txt

set -e -o pipefail

echo "PATH:              ${PATH}"
echo "INSTALL_PATH:      ${INSTALL_PATH}"
echo "VIRAL_NGS_PATH:    ${VIRAL_NGS_PATH}"
echo "MINICONDA_PATH:    ${MINICONDA_PATH}"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV}"

CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"

mkdir -p $INSTALL_PATH/viral-ngs-etc
if [ ! -f $INSTALL_PATH/viral-ngs-etc/viral-ngs ]; then
	ln -s $VIRAL_NGS_PATH $INSTALL_PATH/viral-ngs-etc/viral-ngs
fi
if [ ! -f $INSTALL_PATH/viral-ngs-etc/conda-env ]; then
	ln -s $CONDA_DEFAULT_ENV $INSTALL_PATH/viral-ngs-etc/conda-env
fi

# setup/install viral-ngs directory tree and conda dependencies
sync

# manually install it ourselves instead of using easy-deploy
if [[ "$1" == "minimal" ]]; then
	# a more minimal set of tools (smaller docker image?)
	conda install -y \
		-q $CONDA_CHANNEL_STRING \
		--file "$VIRAL_NGS_PATH/requirements-minimal.txt"
else
	conda install -y \
		-q $CONDA_CHANNEL_STRING \
		--file "$VIRAL_NGS_PATH/requirements-py3.txt" \
		--file "$VIRAL_NGS_PATH/requirements-conda.txt" \
		--file "$VIRAL_NGS_PATH/requirements-conda-tests.txt"
fi

# clean up
conda clean -y --all
