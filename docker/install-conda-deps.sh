#!/bin/bash
#
# Install conda dependencies from one or more requirements files.
#
# Usage: install-conda-deps.sh [options] file1.txt [--x86-only:file2.txt] [file3.txt ...]
#
# Arguments:
#   file.txt              Regular requirements file (installed on all architectures)
#   --x86-only:file.txt   x86-only requirements file (skipped on ARM64)
#
# This script installs all dependencies in a SINGLE resolver call to prevent
# version regressions when derivative images add packages. x86-only files are
# included in the same resolver call on x86, ensuring consistent dependency
# resolution across all packages.
#
# Examples:
#   # Install core packages only
#   install-conda-deps.sh baseimage.txt core.txt
#
#   # Install with x86-only packages (single resolver call on x86, skips x86-only on ARM)
#   install-conda-deps.sh baseimage.txt core.txt --x86-only:core-x86.txt
#
#   # Multiple x86-only files
#   install-conda-deps.sh baseimage.txt core.txt classify.txt --x86-only:classify-x86.txt phylo.txt --x86-only:phylo-x86.txt

set -e -o pipefail

# Architecture detection
ARCH=$(uname -m)

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

# Parse arguments and build requirements list
REQUIREMENTS=""
SKIPPED_X86_FILES=""

for arg in "$@"; do
    if [[ "$arg" == --x86-only:* ]]; then
        # Extract filename from --x86-only:filename.txt
        file="${arg#--x86-only:}"
        if is_x86; then
            if [ -f "$file" ]; then
                echo "Adding x86-only requirements from: $file"
                REQUIREMENTS="$REQUIREMENTS --file $file"
            else
                echo "Warning: x86-only requirements file not found: $file" >&2
            fi
        else
            echo "Skipping x86-only file on $ARCH: $file"
            SKIPPED_X86_FILES="$SKIPPED_X86_FILES $file"
        fi
    else
        # Regular requirements file
        if [ -f "$arg" ]; then
            echo "Adding requirements from: $arg"
            REQUIREMENTS="$REQUIREMENTS --file $arg"
        else
            echo "Warning: requirements file not found: $arg" >&2
        fi
    fi
done

if [ -z "$REQUIREMENTS" ]; then
    if [ -n "$SKIPPED_X86_FILES" ]; then
        echo "No packages to install (all files were x86-only and skipped on $ARCH)"
        exit 0
    fi
    echo "Error: No valid requirements files provided" >&2
    exit 1
fi

echo ""
echo "Installing conda dependencies..."
echo "Architecture: $ARCH"
echo "micromamba version: $(micromamba --version)"
echo ""

# Install all dependencies together in single resolver call
start_keepalive
micromamba install -y -n base $REQUIREMENTS
stop_keepalive

echo ""
echo "Installed packages:"
micromamba list

# Clean up: micromamba clean removes tarballs and index cache, but leaves
# extracted packages in /opt/conda/pkgs/ for hardlink-installed files.
# Trivy scans the full filesystem and flags CVEs in cached JARs/binaries,
# so we nuke the entire pkgs cache to avoid false positives.
micromamba clean -y --all
rm -rf /opt/conda/pkgs/*

# Remove pip's vendoring manifest. pip records the versions it vendored in
# _vendor/vendor.txt, and Trivy's python-pkg analyzer reads that manifest as if each
# line were an installed package -- so `setuptools==70.3.0` in it is reported as
# CVE-2025-47273 even though the env's real setuptools is far newer and the vulnerable
# code (setuptools' PackageIndex download/extract) is not part of what pip vendors.
# The file is a build-time record only; pip runs fine without it (verified: --version,
# list, download, and pkg_resources import all work). Done here rather than per-Dockerfile
# so every tier is covered, including any conda-triggered pip reinstall.
# A find sweep rather than a fixed path so it also catches a repopulated pkgs cache.
find /opt/conda -path '*/pip/_vendor/vendor.txt' -delete 2>/dev/null || true

echo ""
echo "Done installing conda dependencies."
