#!/bin/bash
#
# Install conda dependencies from one or more requirements files.
#
# Usage: install-conda-deps.sh requirements1.txt [requirements2.txt ...]
#
# This script installs all dependencies in a SINGLE resolver call to prevent
# version regressions when derivative images add packages. Always pass ALL
# requirements files together (e.g., core.txt + classify.txt for classify image).

set -e -o pipefail

# Keepalive function for long-running installs (prevents CI timeout)
KEEPALIVE_PID=""

start_keepalive() {
    if [ -n "$KEEPALIVE_PID" ]; then
        return
    fi
    (
        local start_time=$(date +%s)
        while true; do
            sleep 60
            local elapsed=$(($(date +%s) - start_time))
            >&2 echo "Still installing... (${elapsed}s elapsed)"
        done
    ) &
    KEEPALIVE_PID=$!
    disown
}

stop_keepalive() {
    if [ -n "$KEEPALIVE_PID" ]; then
        kill $KEEPALIVE_PID 2>/dev/null || true
        unset KEEPALIVE_PID
    fi
}

trap stop_keepalive EXIT SIGINT SIGQUIT SIGTERM

# Build requirements arguments
REQUIREMENTS=""
for reqfile in "$@"; do
    if [ -f "$reqfile" ]; then
        echo "Adding requirements from: $reqfile"
        REQUIREMENTS="$REQUIREMENTS --file $reqfile"
    else
        echo "Warning: requirements file not found: $reqfile" >&2
    fi
done

if [ -z "$REQUIREMENTS" ]; then
    echo "Error: No valid requirements files provided" >&2
    exit 1
fi

echo ""
echo "Installing conda dependencies..."
echo "micromamba version: $(micromamba --version)"
echo ""

# Install all dependencies together
start_keepalive
micromamba install -y -n base $REQUIREMENTS
stop_keepalive

echo ""
echo "Installed packages:"
micromamba list

# Clean up
micromamba clean -y --all

echo ""
echo "Done installing conda dependencies."
