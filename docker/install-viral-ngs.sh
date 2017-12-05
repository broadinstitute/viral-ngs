#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/viral-ngs)
# and VIRAL_NGS_PATH (typically /opt/viral-ngs/source) to be set.
#
# A miniconda install must exist at /opt/miniconda
#
# Otherwise, this only requires the existence of the following files:
#    easy-deploy-script/easy-deploy-viral-ngs.sh
#    requirements-conda.txt
#    requirements-conda-tests.txt
#    requirements-py3.txt

set -e -o pipefail

export VIRAL_CONDA_ENV_PATH=/opt/miniconda

mkdir -p $INSTALL_PATH/viral-ngs-etc
#mkdir -p $VIRAL_NGS_PATH/.git  # this is needed to make the setup script know we have/will have a git checkout
ln -s $VIRAL_NGS_PATH $INSTALL_PATH/viral-ngs-etc/viral-ngs
ln -s $VIRAL_CONDA_ENV_PATH $INSTALL_PATH/viral-ngs-etc/conda-env
#ln $VIRAL_NGS_PATH/easy-deploy-script/easy-deploy-viral-ngs.sh $INSTALL_PATH

# setup/install viral-ngs directory tree and conda dependencies
sync
#$INSTALL_PATH/easy-deploy-viral-ngs.sh setup-git-local

# manually install it ourselves instead of using easy-deploy
if [[ "$1" == "minimal" ]]; then
	# a more minimal set of tools (smaller docker image?)
	conda install --override-channels -y \
		-q -c broad-viral -c bioconda -c conda-forge -c defaults -c r \
		--file "$VIRAL_NGS_PATH/requirements-minimal.txt"
else
	conda install --override-channels -y \
		-q -c broad-viral -c bioconda -c conda-forge -c defaults -c r \
		--file "$VIRAL_NGS_PATH/requirements-py3.txt" \
		--file "$VIRAL_NGS_PATH/requirements-conda.txt" \
		--file "$VIRAL_NGS_PATH/requirements-conda-tests.txt"
fi

# clean up
conda clean -y --all
