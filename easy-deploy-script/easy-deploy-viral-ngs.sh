#!/bin/bash

#set -e -o pipefail

STARTING_DIR=$(pwd)

# way to get the absolute path to this script that should
# work regardless of whether or not this script has been sourced
# Find original directory of bash script, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        SOURCE="$(readlink "$SOURCE")"
    else
        SOURCE="$(readlink -f "$SOURCE")"
    fi
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT=$SOURCE
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$(echo $SCRIPT_DIRNAME)" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"

CONDA_PREFIX_LENGTH_LIMIT=80

CONTAINING_DIR="viral-ngs-etc"
VIRAL_NGS_DIR="viral-ngs"
CONDA_ENV_BASENAME="conda-env"
CONDA_ENV_CACHE="conda-cache"
PROJECTS_DIR="projects"
MINICONDA_DIR="mc3"

VIRAL_CONDA_ENV_PATH="$SCRIPTPATH/$CONTAINING_DIR/$CONDA_ENV_BASENAME"
PROJECTS_PATH="$SCRIPTPATH/$CONTAINING_DIR/$PROJECTS_DIR"
VIRAL_NGS_PATH="$SCRIPTPATH/$CONTAINING_DIR/$VIRAL_NGS_DIR"
MINICONDA_PATH="$SCRIPTPATH/$CONTAINING_DIR/$MINICONDA_DIR"


# determine if this script has been sourced
# via: http://stackoverflow.com/a/28776166/2328433
([[ -n $ZSH_EVAL_CONTEXT && $ZSH_EVAL_CONTEXT =~ :file$ ]] ||
 [[ -n $KSH_VERSION && $(cd "$(dirname -- "$0")" &&
    printf '%s' "${PWD%/}/")$(basename -- "$0") != "${.sh.file}" ]] ||
 [[ -n $BASH_VERSION && $0 != "$BASH_SOURCE" ]]) && sourced=1 || sourced=0

current_prefix_length=$(echo $MINICONDA_PATH | wc -c | sed -n '1h;1!H;${;g;s/^[ \t]*//g;s/[ \t]*$//g;p;}') # sed trims whitespace
if [ $current_prefix_length -ge $CONDA_PREFIX_LENGTH_LIMIT ]; then
    echo "ERROR: The conda path to be created by this script is too long to work with conda ($current_prefix_length characters):"
    echo "$MINICONDA_PATH"
    echo "This is a known bug in conda ($CONDA_PREFIX_LENGTH_LIMIT character limit): "
    echo "https://github.com/conda/conda-build/pull/877"
    echo "To prevent this error, move this script higher in the filesystem hierarchy."
    exit 1
fi

python_check=$(hash python)
if [ $? -ne 0 ]; then
    echo "It looks like Python is not installed. Exiting."
    if [[ $sourced -eq 0 ]]; then
        exit 1
    else
        return 1
    fi
fi

ram_check=$(python -c "bytearray(768000000)" &> /dev/null)
if [ $? -ne 0 ]; then
    echo ""
    echo "Unable to allocate 768MB."
    echo "=============================================================="
    echo "It appears your current system does not have enough free RAM."
    echo "Consider logging in to a machine with more available memory."
    echo "=============================================================="
    echo ""

    if [[ $sourced -eq 0 ]]; then
        exit 1
    else
        return 1
    fi
fi

function set_locale(){
    export LANG="$1"
    export LC_CTYPE="$1"
    export LC_NUMERIC="$1"
    export LC_TIME="$1"
    export LC_COLLATE="$1"
    export LC_MONETARY="$1"
    export LC_MESSAGES="$1"
    export LC_PAPER="$1"
    export LC_NAME="$1"
    export LC_ADDRESS="$1"
    export LC_TELEPHONE="$1"
    export LC_MEASUREMENT="$1"
    export LC_IDENTIFICATION="$1"
    export LC_ALL="$1"
}

if [[ "$OSTYPE" == "darwin"* ]]; then
    set_locale "en_US.UTF-8"
else
    set_locale "en_US.utf8"
fi

function prepend_miniconda(){
    if [ -d "$MINICONDA_PATH/bin" ]; then
        echo "Miniconda installed."

        echo "Prepending miniconda to PATH..."
        export PATH="$MINICONDA_PATH/bin:$PATH"
        hash -r

        # update to the latest conda this way, since the shell script 
        # is often months out of date
        #if [ -z "$SKIP_CONDA_UPDATE" ]; then
        #    echo "Updating conda..."
        #    conda update -y conda
        #fi
    else
        echo "Miniconda directory not found."
        exit 1
    fi
}

function install_miniconda(){
    if [ -d "$MINICONDA_PATH/bin" ]; then
        echo "Miniconda directory exists."
    else
        echo "Downloading and installing Miniconda..."

        if [[ "$(python -c 'import os; print(os.uname()[0])')" == "Darwin" ]]; then
            miniconda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
        else
            miniconda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        fi

        wget $miniconda_url -O Miniconda3-latest-x86_64.sh -P $(dirname $MINICONDA_PATH)/

        chmod +x $(dirname $MINICONDA_PATH)/Miniconda3-latest-x86_64.sh
        $(dirname $MINICONDA_PATH)/Miniconda3-latest-x86_64.sh -b -f -p "$MINICONDA_PATH"

        rm $(dirname $MINICONDA_PATH)/Miniconda3-latest-x86_64.sh
    fi

    if [ -d "$MINICONDA_PATH/bin" ]; then
        prepend_miniconda
        conda install -q -y conda==4.0.10
        conda install -q -y conda-build==1.7.1
    else
        echo "It looks like the Miniconda installation failed"
        exit 1
    fi
}

function create_project(){
    echo "Populating project directory..."
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
    touch samples-metagenomics.txt
    touch samples-depletion.txt
    touch samples-assembly.txt
    touch samples-runs.txt
    touch samples-assembly-failures.txt
    cp $VIRAL_NGS_PATH/pipes/config.yaml ../../$VIRAL_NGS_DIR/pipes/Snakefile ./
    ln -s "$VIRAL_NGS_PATH/" "$(pwd)/bin"
    ln -s "$VIRAL_CONDA_ENV_PATH/" "$(pwd)/conda-env"
    ln -s "$MINICONDA_PATH/" "$(pwd)/mc3"
    ln -s "$VIRAL_NGS_PATH/pipes/Broad_UGER/run-pipe.sh" "$(pwd)/run-pipe_UGER.sh"
    ln -s "$VIRAL_NGS_PATH/pipes/Broad_LSF/run-pipe.sh" "$(pwd)/run-pipe_LSF.sh"

    cd $starting_dir
}

function activate_env(){
    if [ -d "$SCRIPTPATH/$CONTAINING_DIR" ]; then
        echo "viral-ngs parent directory found"
    else
        echo "viral-ngs parent directory not found: $SCRIPTPATH/$CONTAINING_DIR not found."
        echo "Have you run the setup?"
        echo "Usage: $0 setup"
        cd $STARTING_DIR
        return 1
    fi

    if [ -d "$VIRAL_CONDA_ENV_PATH" ]; then
        if [ -z "$CONDA_DEFAULT_ENV" ]; then
            echo "Activating viral-ngs environment..."
            prepend_miniconda
            
            source activate $VIRAL_CONDA_ENV_PATH

            # unset JAVA_HOME if set, so we can use the conda-supplied version
            if [ ! -z "$JAVA_HOME" ]; then
                unset JAVA_HOME
            fi

            # override $PS1 to have a shorter prompt
            export PS1="(\033[1mviral-ngs\033[0m)\s:\h:\w \! \$ "
        else
            if [[ "$CONDA_DEFAULT_ENV" != "$VIRAL_CONDA_ENV_PATH" ]]; then
                echo "It looks like a conda environment is already active,"
                echo "however it is not the viral-ngs environment."
                echo "To use viral-ngs with your project, deactivate the"
                echo "current environment and then source this file."
                echo "Example: source deactivate && source $(basename $SCRIPT) load"
            else
                echo "The viral-ngs environment is already active."
            fi
            return 0
        fi
    else
        echo "$VIRAL_CONDA_ENV_PATH/ does not exist. Exiting."
        cd $STARTING_DIR
        return 1
    fi
}

function print_usage(){
    echo "Usage: $(basename $SCRIPT) {load,create-project,setup|setup-py2,upgrade}"
}

function symlink_viral_ngs(){
    # the conda-installed viral-ngs folder resides within the
    # opt/ directory of the conda environment, but it contains
    # a version number, so we'll ls and grep for the name
    # and assume it's the first one to show up
    if [ ! -L "$VIRAL_NGS_PATH" ]; then
        echo "Linking to current viral-ngs install..."
        EXPECTED_VIRAL_NGS_VERSION=$(conda list viral-ngs | grep viral-ngs | grep -v packages | awk -F" " '{print $2}')
        VIRAL_NGS_CONDA_PATH="$VIRAL_CONDA_ENV_PATH/opt/"$(ls -1 "$VIRAL_CONDA_ENV_PATH/opt/" | grep "$EXPECTED_VIRAL_NGS_VERSION" | grep -m 1 "viral-ngs")

        if [ -d "$VIRAL_NGS_CONDA_PATH" ]; then
            ln -s "$VIRAL_NGS_CONDA_PATH" "$VIRAL_NGS_PATH"
        else
            echo "Could not find viral-ngs install in conda env:"
            echo "$VIRAL_NGS_CONDA_PATH"
            exit 1
        fi
    else
        echo "$VIRAL_NGS_DIR/ symlink already exists. Skipping link."
    fi
}

function check_viral_ngs_version(){
    # this must be run within an active conda environment
    # so, after a call to "activate_env"
    if [ -z "$SKIP_VERSION_CHECK" ]; then
        echo "Checking viral-ngs version..."
        CURRENT_VER="$(conda list --no-pip viral-ngs | grep viral-ngs | grep -v packages | awk -F" " '{print $2}')"
        # perhaps a better way...
        AVAILABLE_VER="$(conda search --override-channels -f -c r -c bioconda -c conda-forge -c defaults --override-channels viral-ngs --json | grep version | tail -n 1 | awk -F" " '{print $2}' | perl -lape 's/"//g')"
        if [ "$CURRENT_VER" != "$AVAILABLE_VER" ]; then
            echo ""
            echo "============================================================================================================"
            echo "Your copy of viral-ngs appears to be outdated ($CURRENT_VER). A newer version is available ($AVAILABLE_VER)."
            echo "Check the release notes and consider upgrading:"
            echo "https://github.com/broadinstitute/viral-ngs/releases"
            echo "To upgrade: $(basename $SCRIPT) upgrade"
            echo "============================================================================================================"
            echo ""
        else
            echo "viral-ngs is up to date ($CURRENT_VER)"
        fi
    fi
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
       "setup"|"setup-py2")
            if [ $# -eq 1 -o $# -eq 2 ]; then
                if [[ $sourced -eq 1 ]]; then
                    echo "ABORTING. $(basename $SCRIPT) must not be sourced during setup"
                    echo "Usage: $(basename $SCRIPT) setup"
                    return 1
                else
                    if [ ! -z "$CONDA_DEFAULT_ENV" ]; then
                        echo "The viral-ngs setup cannot be run while a conda environment is active."
                        echo "The current environment must first be disabled via 'source deactivate'"
                        exit 1
                    fi

                    mkdir -p $SCRIPTPATH/$CONTAINING_DIR
                    cd $SCRIPTPATH/$CONTAINING_DIR

                    install_miniconda

                    if [ ! -d "$VIRAL_CONDA_ENV_PATH" ]; then
                        # provide an option to use Python 2 in the conda environment
                        if [ "$1" == "setup-py2" ]; then
                            conda create -c r -c bioconda -c conda-forge -c defaults --override-channels -y -p $VIRAL_CONDA_ENV_PATH python=2
                        else
                            conda create -c r -c bioconda -c conda-forge -c defaults --override-channels -y -p $VIRAL_CONDA_ENV_PATH python=3
                        fi
                        
                        # provide an avenue to specify a package path, or to use a previously-built local package
                        if [ $# -eq 2 ]; then
                            if [ "$2" == "--use-local" ]; then
                                conda install -c r -c bioconda -c conda-forge -c defaults --override-channels -y -p $VIRAL_CONDA_ENV_PATH --use-local viral-ngs
                                echo "using local...."
                                exit 0
                            else
                                conda install -c r -c bioconda -c conda-forge -c defaults --override-channels -y -p $VIRAL_CONDA_ENV_PATH $2
                            fi
                        elif [ $# -eq 1 ]; then
                            conda install -c r -c bioconda -c conda-forge -c defaults --override-channels -y -p $VIRAL_CONDA_ENV_PATH viral-ngs
                        fi

                    else
                        echo "$VIRAL_CONDA_ENV_PATH/ already exists. Skipping conda env setup."
                    fi

                    activate_env

                    symlink_viral_ngs

                    # install tools
                    py.test $VIRAL_NGS_PATH/test/unit/test_tools.py

                    # get the version of gatk expected based on the installed conda package
                    EXPECTED_GATK_VERSION=$(conda list | grep gatk | awk -F" " '{print $2}')
                    if [ -z "$GATK_JAR_PATH" ]; then
                        # if the env var is not set, try to get the jar location using the default Broad path
                        if [[ "$(hash dnsdomainname &> /dev/null && dnsdomainname || echo '')" == *"broadinstitute.org" || "$HOSTNAME" == *".broadinstitute.org" || "$DOMAINNAME" == "broadinstitute.org" ]]; then
                            echo "This script is being run on a Broad Institute system."
                            echo "Trying to find GATK..."
                            export GATK_JAR_PATH=$(ls /humgen/gsa-hpprojects/GATK/bin &> /dev/null && sleep 5 && find /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-$EXPECTED_GATK_VERSION-* -maxdepth 0 -type d)/GenomeAnalysisTK.jar
                        fi
                    fi

                    # if the gatk jar file exists, call gatk-register
                    if [ -e "$GATK_JAR_PATH" ]; then
                        echo "GATK found: $GATK_JAR_PATH"
                        gatk-register $GATK_JAR_PATH
                    else
                        echo "GATK jar could not be found on this system for GATK version $EXPECTED_GATK_VERSION"
                        echo "Please activate the viral-ngs conda environment and 'gatk-register /path/to/GenomeAnalysisTK.jar'"
                        exit 0
                    fi

                    echo ""
                    if [ ! -z "$NOVOALIGN_PATH" ]; then
                        novoalign-license-register "$NOVOALIGN_PATH/novoalign.lic"
                    elif [ ! -z "$NOVOALIGN_LICENSE_PATH" ]; then
                        novoalign-license-register "$NOVOALIGN_LICENSE_PATH"
                    else
                        echo "No Novoalign license found via NOVOALIGN_PATH or NOVOALIGN_LICENSE_PATH"
                        echo "Please activate the viral-ngs conda environment and run 'novoalign-license-register /path/to/novoalign.lic'"
                    fi

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
                    activate_env

                    check_viral_ngs_version

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
       "upgrade")
            if [ $# -eq 1 ]; then
                if [[ $sourced -eq 1 ]]; then
                    echo "ABORTING. $(basename $SCRIPT) must not be sourced during upgrade"
                    echo "Usage: $(basename $SCRIPT) upgrade"
                    return 1
                else
                    echo "Upgrading viral-ngs..."

                    if [ -z "$CONDA_DEFAULT_ENV" ]; then
                        activate_env
                    else
                        if [ "$CONDA_DEFAULT_ENV" != "$VIRAL_CONDA_ENV_PATH" ]; then
                            echo "A viral-ngs upgrade cannot be run while a conda environment is active."
                            echo "The current environment must first be disabled via 'source deactivate'"
                            exit 1
                        fi
                    fi

                    if [ ! -z "$(conda list viral-ngs | grep viral-ngs | grep -v packages | awk -F" " '{print $2}')" ]; then
                        if [ -L "$VIRAL_NGS_PATH" ]; then
                            rm $VIRAL_NGS_PATH # remove symlink
                        fi
                        conda update -y -c r -c bioconda -c conda-forge -c defaults --override-channels viral-ngs

                        # recreate symlink to folder for latest viral-ngs in conda-env/opt/
                        symlink_viral_ngs

                        echo ""
                        echo "=================================================================================="
                        echo "Note that if viral-ngs-derived files are present in your project folders,"
                        echo "they may need to be updated. Check the release notes and commit log for more info:"
                        echo "https://github.com/broadinstitute/viral-ngs/releases"
                        echo "=================================================================================="
                        echo ""
                    else
                        echo "The viral-ngs package does not appear to be installed. Have you run setup?"
                        echo "Usage: $(basename $SCRIPT) setup"
                        exit 1
                    fi
                    #source deactivate
                    exit 0
                fi
            else
                echo "Usage: source $(basename $SCRIPT) upgrade"
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

                if [[ "$CONDA_DEFAULT_ENV" != "$VIRAL_CONDA_ENV_PATH" ]]; then
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
