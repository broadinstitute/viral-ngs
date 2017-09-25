#!/bin/sh
set -e -x -o pipefail

pytest test/unit test/integration
