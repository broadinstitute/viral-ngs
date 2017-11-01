#!/bin/bash

set -e

conda install -q -y conda-build==3.0.25
conda install -q -y anaconda-client
conda create -m -q -y -p tools/conda-tools/default Jinja2==2.8
