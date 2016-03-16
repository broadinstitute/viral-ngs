#!/bin/bash
set -e

# the miniconda directory may exist if it has been restored from cache
if [ -d "$MINICONDA_DIR" ] && [ -e "$MINICONDA_DIR/bin/conda" ]; then
    echo "Miniconda install already present from cache: $MINICONDA_DIR"
else # if it does not exist, we need to install miniconda
    rm -rf "$MINICONDA_DIR" # remove the directory in case we have an empty cached directory
    
    if [[ "$TRAVIS_PYTHON_VERSION" == 2* ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi

    bash miniconda.sh -b -p "$MINICONDA_DIR"
    chown -R "$USER" "$MINICONDA_DIR"
    export PATH="$PATH:$MINICONDA_DIR/bin"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda config --add channels bioconda
    conda config --add channels r
    conda update -q conda
    conda info -a # for debugging
fi