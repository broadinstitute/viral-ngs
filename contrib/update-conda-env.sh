#!/bin/bash

function cmd_exists() {
    if hash $@ 2> /dev/null ; then
        return 0
    else
        return 1
    fi
}

# if conda is available
if cmd_exists "conda"; then
    # if a conda environment is active
    if [ ! -z "$CONDA_DEFAULT_ENV" ]; then
        # update the conda environment to install the requirements specified
        conda install -y -c broad-viral -c bioconda --file requirements-conda.txt
    fi
fi
