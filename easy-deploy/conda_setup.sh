#!/bin/bash

set -e -o pipefail

hash -r
conda config --set always_yes yes --set changeps1 no
conda config --add channels bioconda
conda config --add channels r
conda update -q conda