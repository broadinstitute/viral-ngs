#!/bin/bash
#
# Install conda dependencies from one or more requirements files.
#
# Usage: install-conda-deps.sh [options] requirements1.txt [requirements2.txt ...]
#
# Options:
#   --x86-only    Only install if running on x86_64 architecture (skip on ARM)
#
# This script installs all dependencies in a SINGLE resolver call to prevent
# version regressions when derivative images add packages. Always pass ALL
# requirements files together (e.g., core.txt + classify.txt for classify image).

set -e -o pipefail

# Architecture detection
ARCH=$(uname -m)
X86_ONLY=false

is_x86() {
    [[ "$ARCH" == "x86_64" || "$ARCH" == "amd64" ]]
}

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

# Parse options and build requirements arguments
REQUIREMENTS=""
for arg in "$@"; do
    case "$arg" in
        --x86-only)
            X86_ONLY=true
            ;;
        *)
            if [ -f "$arg" ]; then
                echo "Adding requirements from: $arg"
                REQUIREMENTS="$REQUIREMENTS --file $arg"
            else
                echo "Warning: requirements file not found: $arg" >&2
            fi
            ;;
    esac
done

# Skip x86-only packages on non-x86 architectures
if $X86_ONLY && ! is_x86; then
    echo "Skipping x86-only packages on $ARCH architecture"
    exit 0
fi

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
