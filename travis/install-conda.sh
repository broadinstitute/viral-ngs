#!/bin/bash
set -e

# the miniconda directory may exist if it has been restored from cache
# if it does not exist, we need to install miniconda
if [ ! -d $MINICONDA_DIR ]; then
    if [[ "$TRAVIS_PYTHON_VERSION" == 2* ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi

    bash miniconda.sh -b -p $MINICONDA_DIR
    chown -R $USER $MINICONDA_DIR
    export PATH="$PATH:$MINICONDA_DIR/bin"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda config --add channels bioconda
    conda config --add channels r
    conda update -q conda
    conda info -a # for debugging
fi