#!/bin/bash

apt-get install -y -qq --no-install-recommends gcc
conda install -q -y conda-build anaconda-client Jinja2
