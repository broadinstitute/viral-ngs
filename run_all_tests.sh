#!/bin/sh
set -e -x -o pipefail

py.test test/unit test/integration
