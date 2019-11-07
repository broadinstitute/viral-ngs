#!/bin/bash

set -e -o pipefail

pushd docs
make html && echo "Docs built successfully!" || echo "Docs did NOT build successfully."
popd
