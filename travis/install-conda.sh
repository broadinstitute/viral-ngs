#!/bin/bash
set -e

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b -p $HOME/miniconda
chown -R $USER $HOME/miniconda
export PATH="$PATH:$HOME/miniconda/bin"
hash -r
conda config --set always_yes yes --set changeps1 no
conda config --add channels bioconda
conda config --add channels r
conda update -q conda
conda info -a # for debugging