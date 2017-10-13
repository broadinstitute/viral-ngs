#!/bin/bash
set -e

# the miniconda directory may exist if it has been restored from cache
if [ -d "$MINICONDA_DIR" ] && [ -e "$MINICONDA_DIR/bin/conda" ]; then
    echo "Miniconda install already present from cache: $MINICONDA_DIR"
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        # on OSX we need to rely on the conda Python rather than the Travis-supplied system Python
        # so conda has a higher precedence
        export PATH="$MINICONDA_DIR/bin:$PATH"
    else
        export PATH="$MINICONDA_DIR/bin:$PATH"
    fi
    hash -r
else # if it does not exist, we need to install miniconda
    rm -rf "$MINICONDA_DIR" # remove the directory in case we have an empty cached directory

    if [[ "$TRAVIS_PYTHON_VERSION" == 2* ]]; then
        if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
            wget --quiet https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh -O miniconda.sh;
        else
            wget --quiet https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
        fi
     else
        if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
            wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
        else
            wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
        fi
    fi

    bash miniconda.sh -b -p "$MINICONDA_DIR"
    chown -R "$USER" "$MINICONDA_DIR"
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        # on OSX we need to rely on the conda Python rather than the Travis-supplied system Python
        # so conda has a higher precedence
        export PATH="$MINICONDA_DIR/bin:$PATH"
    else
        export PATH="$MINICONDA_DIR/bin:$PATH"
    fi
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda config --set anaconda_upload yes # for uploading packages after successful build
    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install --quiet -y conda #conda=4.2 # pin to 4.2.* until this is fixed: https://github.com/conda/conda-build/issues/1666
    conda config --set auto_update_conda false
    conda install --quiet -y java-jdk==8.0.112
    conda install --quiet -y conda-build # needed to build recipe
    conda install --quiet -y anaconda-client # needed to upload build package to anaconda.org
    #conda install --quiet -y -c conda-forge -f curl # the bioconda curl is broken as of 21 Feb 2017
fi

conda info -a # for debugging
