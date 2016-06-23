#!/bin/bash

STARTING_DIR=$(pwd)

# way to get the absolute path to this script that should
# work regardless of whether or not this script has been sourced
SCRIPT="$(readlink --canonicalize-existing "${BASH_SOURCE[0]}")"
SCRIPTPATH="$(dirname "$SCRIPT")"

CONTAINING_DIR="viral-ngs-analysis-software"
VIRAL_NGS_DIR="viral-ngs"
CONDA_ENV_BASENAME="conda-env"
CONDA_ENV_CACHE="conda-cache"
PROJECTS_DIR="projects"

VIRAL_CONDA_ENV_PATH="$SCRIPTPATH/$CONTAINING_DIR/$CONDA_ENV_BASENAME"
VIRAL_CONDA_CACHE_PATH="/broad/hptmp/$(whoami)/$CONTAINING_DIR/$CONDA_ENV_CACHE"
PROJECTS_PATH="$SCRIPTPATH/$CONTAINING_DIR/$PROJECTS_DIR"
VIRAL_NGS_PATH="$SCRIPTPATH/$CONTAINING_DIR/$VIRAL_NGS_DIR"

# determine if this script has been sourced
# via: http://stackoverflow.com/a/28776166/2328433
([[ -n $ZSH_EVAL_CONTEXT && $ZSH_EVAL_CONTEXT =~ :file$ ]] ||
 [[ -n $KSH_VERSION && $(cd "$(dirname -- "$0")" &&
    printf '%s' "${PWD%/}/")$(basename -- "$0") != "${.sh.file}" ]] ||
 [[ -n $BASH_VERSION && $0 != "$BASH_SOURCE" ]]) && sourced=1 || sourced=0

function load_dotkits(){
    source /broad/software/scripts/useuse
    reuse .anaconda3-4.0.0
    
    if [ -z "$GATK_PATH" ]; then
        reuse .gatk3-2.2
        # the Broad sets an alias for GATK, so we need to parse out the path
        export GATK_PATH="$(dirname $(alias | grep GenomeAnalysisTK | perl -lape 's/(.*)\ (\/.*.jar).*/$2/g'))"
        export GATK_JAR="$(alias | grep GenomeAnalysisTK | perl -lape 's/(.*)\ (\/.*.jar).*/$2/g')"
    else
        echo "GATK_PATH is set to '$GATK_PATH'"
        echo "Continuing..."
    fi

    if [ -z "$NOVOALIGN_PATH" ]; then
        reuse .novocraft-3.02.08
        export NOVOALIGN_PATH="$(dirname $(which novoalign))"
    else
        echo "NOVOALIGN_PATH is set to '$NOVOALIGN_PATH'"
        echo "Continuing..."
    fi
}

function create_project(){
    # first arg is project folder name
    starting_dir=$(pwd)

    mkdir -p $PROJECTS_PATH
    cd $PROJECTS_PATH
    mkdir $1
    cd $1
    mkdir data log reports tmp
    cd data
    mkdir 00_raw 01_cleaned 01_per_sample 02_align_to_self 02_assembly 03_align_to_ref 03_interhost 04_intrahost
    cd ../
    touch samples-depletion.txt
    touch samples-assembly.txt
    touch samples-runs.txt
    touch samples-assembly-failures.txt
    cp $VIRAL_NGS_PATH/pipes/config.yaml ../../$VIRAL_NGS_DIR/pipes/Snakefile ./
    ln -s $VIRAL_NGS_PATH/ $(pwd)/bin
    ln -s $VIRAL_CONDA_ENV_PATH/ $(pwd)/venv
    ln -s $VIRAL_NGS_PATH/pipes/Broad_UGER/run-pipe.sh $(pwd)/run-pipe_UGER.sh
    ln -s $VIRAL_NGS_PATH/pipes/Broad_LSF/run-pipe.sh $(pwd)/run-pipe_LSF.sh

    cd $starting_dir
}

function activate_env(){
    if [ -d "$VIRAL_CONDA_ENV_PATH" ]; then
        source activate $VIRAL_CONDA_ENV_PATH
    else
        echo "$VIRAL_CONDA_ENV_PATH/ does not exist. Exiting."
        cd $STARTING_DIR
        return 1
    fi
}

function activate_environment(){
    load_dotkits


    # create the conda cache path if it does not exist
    if [ ! -d "$VIRAL_CONDA_ENV_PATH" ]; then
        mkdir -p $VIRAL_CONDA_CACHE_PATH
    fi
    export CONDA_ENVS_PATH="$VIRAL_CONDA_CACHE_PATH"

    echo "$SCRIPTPATH/$CONTAINING_DIR"
    if [ -d "$SCRIPTPATH/$CONTAINING_DIR" ]; then
        cd $SCRIPTPATH/$CONTAINING_DIR
    else
        echo "viral-ngs parent directory not found: $CONTAINING_DIR not found."
        echo "Have you run the setup?"
        echo "Usage: $0 setup"
        cd $STARTING_DIR
        return 1
    fi

    activate_env
}

function print_usage(){
    echo "Usage: $(basename $SCRIPT) {load,create-project,setup}"
}

if [ $# -eq 0 ]; then
    print_usage
    if [[ $sourced -eq 0 ]]; then
        exit 1
    else
        return 1
    fi
else
    case "$1" in
       "setup")
            if [ $# -eq 1 ]; then
                if [[ $sourced -eq 1 ]]; then
                    echo "ABORTING. $(basename $SCRIPT) must not be sourced during setup"
                    echo "Usage: $(basename $SCRIPT) setup"
                    return 1
                else
                    mkdir -p $SCRIPTPATH/$CONTAINING_DIR
                    cd $SCRIPTPATH/$CONTAINING_DIR

                    # clone viral-ngs if it does not already exist
                    # TODO: avoid duplication; symlink to conda version in opt/
                    #  ...follow symlink to illumina.py after install
                    #  of conda package and link dirname to $VIRAL_NGS_DIR
                    # if symlink does not already exist
                    if [ ! -d "$VIRAL_NGS_PATH" ]; then
                        git clone --depth 5 https://github.com/broadinstitute/viral-ngs.git
                    else
                        echo "$VIRAL_NGS_DIR/ already exists. Skipping clone."
                    fi

                    load_dotkits

                    export CONDA_ENVS_PATH="$VIRAL_CONDA_CACHE_PATH"
                    conda config --add channels bioconda

                    if [ ! -d "$VIRAL_CONDA_ENV_PATH" ]; then
                        conda create -y -p $VIRAL_CONDA_ENV_PATH viral-ngs
                        PYVER=`python -V 2>&1 | cut -c 8`
                        if [ "$PYVER" = "3" ]; then
                            conda remove futures # uninstall futures if on Python3 since sometimes it sneaks in
                        fi
                    else
                        echo "$VIRAL_CONDA_ENV_PATH/ already exists. Skipping python venv setup."
                    fi

                    activate_env

                    #pip install -r $VIRAL_NGS_PATH/requirements.txt
                    #pip install -r $VIRAL_NGS_PATH/requirements-tests.txt
                    #pip install -r $VIRAL_NGS_PATH/requirements-pipes.txt

                    # install tools
                    #py.test $VIRAL_NGS_PATH/test/unit/test_tools.py
                    $VIRAL_NGS_PATH/install_tools.py
                    # TODO: add older versions of GATK to bioconda
                    # or find newer GATK on the Broad cluster
                    gatk-register $GATK_JAR

                    echo "Setup complete. Do you want to start a project? Run:"
                    echo "$0 create-project <project_name>"
                    echo ""
                fi
            else
                echo "Usage: $(basename $SCRIPT) setup"
                if [[ $sourced -eq 0 ]]; then
                    exit 1
                else
                    return 1
                fi
            fi
       ;;
       "load")
            if [ $# -eq 1 ]; then
                if [[ $sourced -eq 0 ]]; then
                    echo "ABORTING. $(basename $SCRIPT) must be sourced."
                    echo "Usage: source $(basename $SCRIPT) load"
                else
                    activate_environment
                    ls -lah
                    return 0
                fi
            else
                echo "Usage: source $(basename $SCRIPT) load"
                if [[ $sourced -eq 0 ]]; then
                    exit 1
                else
                    return 1
                fi
            fi
       ;;
       "create-project")
            if [ $# -ne 2 ]; then
                echo "Usage: $(basename $SCRIPT) create-project <project_name>"
                if [[ $sourced -eq 0 ]]; then
                    exit 1
                else
                    return 1
                fi
            else
                if [ ! -d "$PROJECTS_PATH/$2" ]; then
                    create_project $2 && echo "Project created: $PROJECTS_PATH/$2" && echo "OK"
                else
                    echo "WARNING: $PROJECTS_PATH/$2/ already exists."
                    echo "Skipping project creation."
                fi

                echo ""

                if [[ "$VIRTUAL_ENV" != "$VIRAL_CONDA_ENV_PATH" ]]; then
                    echo "It looks like the vial-ngs environment is not active."
                    echo "To use viral-ngs with your project, source this file."
                    echo "Example: source $(basename $SCRIPT) load"
                else
                    # if the viral-ngs environment is active and we have sourced this file
                    if [[ $sourced -eq 1 ]]; then
                        # change to the project directory
                        if [ -d "$PROJECTS_PATH/$2" ]; then
                            cd "$PROJECTS_PATH/$2"
                        fi
                        return 0
                    fi
                fi
            fi
       ;;
       *)
            print_usage
        ;;
    esac
fi
