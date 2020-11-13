#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/viral-ngs),
# VIRAL_NGS_PATH (typically /opt/viral-ngs/source), and
# CONDA_DEFAULT_ENV (typically /opt/miniconda) to be set.
#
# A miniconda install must exist at $CONDA_DEFAULT_ENV
# and $CONDA_DEFAULT_ENV/bin must be in the PATH

set -e -o pipefail

CONDA_INSTALL_TIMEOUT="90m"

echo "PATH:              ${PATH}"
echo "INSTALL_PATH:      ${INSTALL_PATH}"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "VIRAL_NGS_PATH:    ${VIRAL_NGS_PATH}"
echo "MINICONDA_PATH:    ${MINICONDA_PATH}"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV}"
CONDA_CHANNEL_STRING="--override-channels -c broad-viral -c conda-forge -c bioconda -c defaults"

# solving the dependency graph for a conda environment can take a while.
# so long, in fact, that the conda process can run for >10 minutes without
# writing to stderr/stdout. That means Travis CI is likely to kill the job
# these functions write out periodically to keep the build job alive
#   adapted from:
#     https://github.com/matthew-brett/multibuild/blob/d1252d15e95712700865fe4a3e7c20f978efba03/common_utils.sh
# similar to travis_wait, but with output
#   see: 
#     https://docs.travis-ci.com/user/common-build-problems/#build-times-out-because-no-output-was-received

# Work round bug in travis xcode image described at
# https://github.com/direnv/direnv/issues/210
shell_session_update() { :; }
unset -f cd
unset -f pushd
unset -f popd

function start_keepalive {
    if [ -n "$KEEPALIVE_PID" ]; then
        return
    fi

    >&2 echo "Running..."
    # Start a process that runs as a keep-alive
    # to avoid travis quitting if there is no output
    (while true; do
        sleep 120
        >&2 echo "Still running..."
    done) &
    KEEPALIVE_PID=$!
    disown
}

function stop_keepalive {
    if [ ! -n "$KEEPALIVE_PID" ]; then
        return
    fi

    kill $KEEPALIVE_PID
    unset KEEPALIVE_PID

    >&2 echo "Done."
}
trap stop_keepalive EXIT SIGINT SIGQUIT SIGTERM

# setup/install viral-ngs directory tree and conda dependencies
sync

REQUIREMENTS=""
for condafile in $*; do
	REQUIREMENTS="$REQUIREMENTS --file $condafile"
done

# run conda install with keepalive subshell process running in background
# to keep travis build going. Enforce a hard timeout via timeout GNU coreutil
start_keepalive
#timeout $CONDA_INSTALL_TIMEOUT conda install -y -q $CONDA_CHANNEL_STRING -p "${CONDA_PREFIX}" $REQUIREMENTS
mamba install -y -q $CONDA_CHANNEL_STRING -p "${CONDA_PREFIX}" $REQUIREMENTS
stop_keepalive

# clean up
conda clean -y --all
