#!/bin/bash

set -e -o pipefail

pushd docs
make html
popd
