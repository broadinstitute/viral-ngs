#!/bin/bash

function cmd_exists() {
    if hash $@ 2> /dev/null ; then
        return 0
    else
        return 1
    fi
}

requirements_file="$(mktemp)"
trap "rm $requirements_file" EXIT

# if conda is available
if cmd_exists "conda"; then
    # if a conda environment is active
    if [ ! -z "$CONDA_DEFAULT_ENV" ]; then
        # update the conda environment to install the requirements specified
        if [ "$(uname)" == "Darwin" ]; then
            # some conda packages dont exist on OSX
            cat requirements-conda.txt | grep -v kraken > "$requirements_file"
        else
            # for linux, just use requirements-conda as-is
            cp requirements-conda.txt "$requirements_file"
        fi

        conda install -y -c bioconda --file "$requirements_file"
    fi
fi
