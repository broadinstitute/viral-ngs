#!/bin/bash
set -e -o pipefail

ln -s "$(which gatk3)" "$(dirname $(which gatk3))/gatk"