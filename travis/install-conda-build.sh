#!/bin/bash

apt-get install -y -qq --no-install-recommends gcc
conda install -q -y conda-build==3.17.8 anaconda-client Jinja2==2.8
