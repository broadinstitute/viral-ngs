#!/bin/bash
# This script primarily enables or disables the Travis dependency
# cache depending on whether we're on the master branch or not.
set -e

# Report how big things are
echo "Docker cache space usage:"
du -hs $CACHE_DIR/*
