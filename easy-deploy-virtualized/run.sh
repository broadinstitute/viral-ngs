#!/bin/bash

# 2015 Christopher Tomkins-Tinch
# Broad Institute of MIT and Harvard

function ask() {
    while true; do
 
        if [ "${2:-}" = "Y" ]; then
            prompt="Y/n"
            default=Y
        elif [ "${2:-}" = "N" ]; then
            prompt="y/N"
            default=N
        else
            prompt="y/n"
            default=
        fi
 
        # Ask the question
        read -p "$1 [$prompt] " REPLY
 
        # Default?
        if [ -z "$REPLY" ]; then
            REPLY=$default
        fi
 
        # Check if the reply is valid
        case "$REPLY" in
            Y*|y*) echo " "; return 0 ;;
            N*|n*) echo " "; return 1 ;;
        esac
 
    done
}

function cmd_exists() {
    if hash $@ 2> /dev/null ; then
        return 0
    else
        return 1
    fi
}

function all_commands_exist {
    # Pass in a space-delimited string of commands
    # to try. Returns 0 if all exist
    #
    # $1 array of commands to try

    all_exist=true

    cmds_array=($(echo $1 | tr " " "\n"))

    for f in "${cmds_array[@]}"; do
        if ! cmd_exists $f; then
            echo "Command missing: $f"
            if $all_exist; then
                all_exist=false
            fi
        fi
    done

    if $all_exist; then
        return 0
    else
        return 1
    fi
}

function install_with_blocking() {
    # Use to run an install command where there isn't an easy way
    # to see if the command already exists
    #
    # $1 the install command to run
    # $2 description of thing to be installed

    if eval $1 ; then
        return 0
    else
        echo "Install of $2 seems to have been unsuccessful"
        exit 1
    fi
}

function use_or_install_with_checks() {
    # Used to run an install command if a specified command does
    # not exist
    #
    # $1 command to check to see if it is already installed
    # $2 description of the package
    # $3 command to install it    
    # $4 unattended if 1

    if ! cmd_exists $1  ; then
        if [[ $4 == 1 ]] || $(ask "Install $2?" Y); then
            echo "Installing $2..."
            echo "$3"
            eval $3
            if cmd_exists $1  ; then
                return 0
            else
                echo "Attempt to install $2 failed. Please check error messages and logs."
                exit 1
            fi
        else
            echo "Exiting. Please install $2 manually."
            exit 1
        fi
    else
        echo "$2 found. Continuing..."
        return 0
    fi
}

# in the event that the script is sent a kill signal
# suspend the vm cleanly
function clean_up {
    if cmd_exists "vagrant"; then
        vagrant suspend
    fi
    exit $1
}
trap clean_up SIGHUP SIGINT SIGTERM

echo $' +-----------------------------------------------------------+ '
echo $' |       _   _ _           _        _   _ _____  _____       | '
echo $' |      | | | (_)         | |      | \ | |  __ \/  ___|      | '
echo $' |      | | | |_ _ __ __ _| |______|  \| | |  \/\ `--.       | '
echo $' |      | | | | | \'__/ _` | |______| . ` | | __  `--. \      | '
echo $' |      \ \_/ / | | | (_| | |      | |\  | |_\ \/\__/ /      | '
echo $' |       \___/|_|_|  \__,_|_|      \_| \_/\____/\____/       | '
echo $' |                                                           | '
echo $' |                                                           | '
echo $' |   _____                  ______           _               | '
echo $' |  |  ___|                 |  _  \         | |              | '
echo $' |  | |__  __ _ ___ _   _   | | | |___ _ __ | | ___  _   _   | '
echo $' |  |  __|/ _` / __| | | |  | | | / _ \ \'_ \| |/ _ \| | | |  | '
echo $' |  | |__| (_| \__ \ |_| |  | |/ /  __/ |_) | | (_) | |_| |  | '
echo $' |  \____/\__,_|___/\__, |  |___/ \___| .__/|_|\___/ \__, |  | '
echo $' |                   __/ |            | |             __/ |  | '
echo $' |                  |___/             |_|            |___/   | '
echo $' |                                                           | '
echo $' +-----------------------------------------------------------+ '
echo  '    Hosts known to work: OSX 10.10 and Ubuntu 15.04 x64        '

# cd to directory of this file so we are in the same location as the Vagrantfile
cd ${0%/*}

relevant_commands='ansible vagrant VirtualBox'

if ! [ $(uname -m) == 'x86_64' ]; then
    echo "This script is intended to be run on 64-bit machines."
    exit 1
fi

if ! all_commands_exist "$relevant_commands"; then
    if $(ask "Install local dependencies for running viral-ngs?" Y); then
        echo "Setting up dependencies..."

        if [ -z "$GATK_PATH" ]; then 
            echo "Prior to setting up everything, you must license and download GATK ( https://www.broadinstitute.org/gatk/download/ )." 
            echo "Once downloaded, extract the archive and set the environment variable GATK_PATH to the directory path"
            echo "containing 'GenomeAnalysisTK.jar'"
            echo "  ex. 'export GATK_PATH=/path/to/gatk'"
            echo "       Note: add the export line to .bashrc or .bash_profile for persistence across sessions."
            exit 1;
        else
            echo "GATK_PATH is set to '$GATK_PATH'"
            echo "Continuing..."
        fi

        if [ -z "$NOVOALIGN_PATH" ]; then 
            echo "Prior to setting up everything, you must license and download NovoAlign ( http://www.novocraft.com/products/novoalign/ )." 
            echo "Once downloaded, extract the archive and set the environment variable NOVOALIGN_PATH to the directory path"
            echo "containing 'novoalign'"
            echo "  ex. 'export NOVOALIGN_PATH=/path/to/novoalign'"
            echo "       Note: add the export line to .bashrc or .bash_profile for persistence across sessions."
            exit 1;
        else
            echo "NOVOALIGN_PATH is set to '$NOVOALIGN_PATH'"
            echo "Continuing..."
        fi

        if [ "$(uname)" == "Darwin" ]; then
            echo "Platform: Mac"

            if use_or_install_with_checks "brew" "Homebrew package mananger" 'ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" ' 1 ; then
                echo "Homebrew found, installing cask..."
                use_or_install_with_checks 'ansible' "Ansible deploy and config manager" 'brew install ansible' 1
                if use_or_install_with_checks 'brew-cask' "Homebrew cask" 'brew install caskroom/cask/brew-cask' 1 ; then
                    use_or_install_with_checks 'VirtualBox' "VirtualBox VM manager" 'brew cask install virtualbox' 1
                    use_or_install_with_checks 'vagrant' "Vagrant VM provisioning tool" 'brew cask install vagrant' 1
                    install_with_blocking 'brew cask install vagrant-manager' 'vagrant-manager'
                fi          
                # if the version of ruby is not 2.0 ( OSX < 10.10 )
                if ! echo $(ruby -e 'print RUBY_VERSION') | grep "2." &> /dev/null ; then
                    echo "Ruby is not at version 2.x. Installing 2.0 from Homebrew..."
                    install_with_blocking "brew install ruby20" "Ruby 2.0"
                fi
            fi

        elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
            echo "Platform: Linux"

            if cmd_exists "apt-get"; then
                use_or_install_with_checks "VirtualBox" "VirtualBox VM manager" 'sudo apt-get -y install virtualbox virtualbox-dkms' 1
                #use_or_install_with_checks "ansible" "Ansible deploy and config manager" 'sudo apt-get -y install ansible' 1
                use_or_install_with_checks "ansible" "ansible" "sudo apt-get -y install software-properties-common && sudo apt-add-repository -y ppa:ansible/ansible && sudo apt-get -y update && sudo apt-get -y install ansible" 1
                #use_or_install_with_checks "vagrant" "Vagrant VM provisioning tool" 'sudo apt-get -y install vagrant' 1
                if use_or_install_with_checks "wget" "wget" "sudo apt-get -y install wget" 1; then                        
                    use_or_install_with_checks "vagrant" "vagrant" "wget https://dl.bintray.com/mitchellh/vagrant/vagrant_1.7.4_x86_64.deb && sudo dpkg --install vagrant_1.7.4_x86_64.deb && rm vagrant_1.7.4_x86_64.deb" 1
                fi
                install_with_blocking 'sudo apt-get -y install zlibc zlib1g zlib1g-dev' 'zlibc'

                if use_or_install_with_checks "ruby" "Ruby" 'sudo apt-get -y install ruby' 1 ; then                
                    if ! echo $(ruby -e 'print RUBY_VERSION') | grep "2." &> /dev/null ; then
                        echo "Please update your system ruby to >2.0 and re-run this script"
                        exit 1
                    fi
                fi
            else
                echo "This script is intended to be used on a system with apt. If you are using a Linux distro without apt, please set up the environment manually."
                exit 1
            fi

        elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
            echo "Platform: Windows. Not currently supported."
            exit 1
        fi
    else
        echo "Exiting."
        exit 1
    fi
fi

if all_commands_exist $relevant_commands; then
    echo $'Dependencies for building the viral-ngs environment are SATISFIED.\n'
    echo $'Next you have an opportunity to choose whether you would like to provision'
    echo $'the viral-ngs environment LOCALLY or on AWS EC2.'
    echo $'In either case, you will be directed to an ssh session on a VM where'
    echo $'you can work with viral-ngs.\n'

    read -p "Run viral-ngs locally or on AWS EC2? [LOCAL/aws] " REPLY
    
    # Default?
    if [ -z "$REPLY" ]; then
        REPLY="local"
    fi
    
    # Check if the reply is valid
    case "$REPLY" in
        L*|l*) # could also be explicit LOCAL|local)
            echo "Running locally..."; 
            vagrant up --provider=virtualbox
            
            if ! [ $? -eq 0 ]; then
                echo "The VM was not set up correctly. Check above for error messages."
                exit 1
            fi

            vagrant box update
            vagrant ssh
            
            vagrant suspend
            if [ $? -eq 0 ]; then
                echo $'\nThe Viral-NGS virtual machine has been suspended. Next time you run this'
                echo $'script, you will be connected exactly where you left off.\n'
                exit 0
            else
                echo $'The Vial-NGS virtual machine did not close cleanly.'
                exit 1
            fi
            ;;
        A*|a*)
            echo "Starting on AWS...";

            if ! echo $(ruby -e 'print RUBY_VERSION') | grep "2." &> /dev/null ; then
                echo "Running viral-ngs on AWS needs ruby >2.0"
                exit 1
            else
                if ! echo $(vagrant plugin list) | grep "aws" &> /dev/null ; then
                    install_with_blocking "vagrant plugin install vagrant-aws" 'vagrant-aws'
                fi
            fi

            if [ -z "$EC2_ACCESS_KEY_ID" ] || [ -z "$EC2_SECRET_ACCESS_KEY" ] || [ -z "$EC2_KEYPAIR_NAME" ] || [ -z "$EC2_PRIVATE_KEY_PATH" ] || [ -z "$EC2_REGION" ] || [ -z "$EC2_SECURITY_GROUP" ]; then 
                echo "==========================================================================================="
                echo "Prior to deploying to EC2 you must obtain AWS IAM credentials permitting instance creation." 
                echo "You must also create a private/public keypair for use with EC2." 
                echo "The following environment variables are not set:"
                echo ""
                if [ -z "$EC2_ACCESS_KEY_ID" ]; then
                    echo "    EC2_ACCESS_KEY_ID"
                fi
                if [ -z "$EC2_SECRET_ACCESS_KEY" ]; then
                    echo "    EC2_SECRET_ACCESS_KEY"
                fi
                if [ -z "$EC2_REGION" ]; then
                    echo "    EC2_REGION"
                fi
                if [ -z "$EC2_KEYPAIR_NAME" ]; then
                    echo "    EC2_KEYPAIR_NAME"
                fi
                if [ -z "$EC2_PRIVATE_KEY_PATH" ]; then
                    echo "    EC2_PRIVATE_KEY_PATH"
                fi
                if [ -z "$EC2_SECURITY_GROUP" ]; then
                    echo "    EC2_SECURITY_GROUP"
                fi
                echo ""
                echo "If you wish, you may specify these values now for use in only this session."
                echo "Unless the environment variables are set you will be prompted again next time"
                echo "you run this script."
                echo ""
                if $(ask "Specify the missing values now?" Y); then
                    while [ -z "$EC2_ACCESS_KEY_ID" ]; do
                        read -p "What is the AWS access_key_id to use? " EC2_ACCESS_KEY_ID
                        if [ $(echo ${#EC2_ACCESS_KEY_ID}) -gt 0 ]; then
                            export EC2_ACCESS_KEY_ID="$EC2_ACCESS_KEY_ID"
                        fi
                    done

                    while [ -z "$EC2_SECRET_ACCESS_KEY" ]; do
                        read -p "What is the AWS secret_access_key to use? " EC2_SECRET_ACCESS_KEY
                        if [ $(echo ${#EC2_SECRET_ACCESS_KEY}) -gt 0 ]; then
                            export EC2_SECRET_ACCESS_KEY="$EC2_SECRET_ACCESS_KEY"
                        fi
                    done

                    while [ -z "$EC2_REGION" ]; do
                        read -p "What is the EC2 region (ex. 'us-east-1')? " EC2_REGION
                        if [ $(echo ${#EC2_REGION}) -gt 0 ]; then
                            export EC2_REGION="$EC2_REGION"
                        fi
                    done

                    while [ -z "$EC2_KEYPAIR_NAME" ]; do
                        read -p "What is AWS keypair name (region-specific)? " EC2_KEYPAIR_NAME
                        if [ $(echo ${#EC2_KEYPAIR_NAME}) -gt 0 ]; then
                            export EC2_KEYPAIR_NAME="$EC2_KEYPAIR_NAME"
                        fi
                    done

                    while [ -z "$EC2_PRIVATE_KEY_PATH" ]; do
                        read -p "Where is the local private key (*.pem) file for the AWS keypair (file path)? " EC2_PRIVATE_KEY_PATH
                        if [ $(echo ${#EC2_PRIVATE_KEY_PATH}) -gt 0 ]; then
                            export EC2_PRIVATE_KEY_PATH="$EC2_PRIVATE_KEY_PATH"
                        fi
                    done

                    while [ -z "$EC2_SECURITY_GROUP" ]; do
                        read -p "Specify the *name* of the EC2 security group to use; the group must permit access on TCP port 22 (region-specific): " EC2_SECURITY_GROUP
                        if [ $(echo ${#EC2_SECURITY_GROUP}) -gt 0 ]; then
                            export EC2_SECURITY_GROUP="$EC2_SECURITY_GROUP"
                        fi
                    done
                else
                    echo "Exiting. Please set the AWS environment variables."
                    exit 1
                fi
            else
               echo "EC2_ACCESS_KEY_ID is set."
               echo "EC2_SECRET_ACCESS_KEY is set."
               echo "EC2_KEYPAIR_NAME is set."
               echo "EC2_PRIVATE_KEY_PATH is set."
               echo "EC2_REGION is set."
               echo "EC2_SECURITY_GROUP is set."
               echo "Continuing..."
            fi

            echo "====================================================================================="
            echo "Data sync to the AWS VM is currently one-way. You will be responsible for downloading"
            echo "or otherwise transferring data from the AWS VM once processing is complete."
            echo ""
            if $(ask "Have you copied any new data you would like to work with in this session to 'data/' ?" Y); then
                vagrant up --provider=aws

                if ! [ $? -eq 0 ]; then
                    echo "The VM was not set up correctly. Check above for error messages."
                    vagrant destroy
                    exit 1
                fi

                vagrant rsync # one-way sync of files, local->ec2
                vagrant ssh
            fi

            ;;
        *)
            echo "Exiting."
        ;;
    esac
fi
