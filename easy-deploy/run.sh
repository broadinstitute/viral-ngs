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

function install_with_blocking() {
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
            echo "Aborting. Please install $2 manually."
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

#echo $'viral-ngs easy deployment script\n'
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


if ! cmd_exists "vagrant" || ! cmd_exists "VirtualBox"; then
    if $(ask "Install local dependencies for running viral-ngs?" Y); then
        echo "Setting up dependencies..."

        if [ "$(uname)" == "Darwin" ]; then
            echo "Platform: Mac"

            if use_or_install_with_checks "brew" "Homebrew package mananger" 'ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" ' 1 ; then
                echo "Homebrew found, installing cask..."
                if use_or_install_with_checks 'brew-cask' "Homebrew cask" 'brew install caskroom/cask/brew-cask' 1 ; then
                    install_with_blocking 'brew cask install virtualbox'
                    install_with_blocking 'brew cask install vagrant'
                    install_with_blocking 'brew cask install vagrant-manager'
                fi        
            fi

        elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
            echo "Platform: Linux"

            if [ $(python -mplatform | grep Ubuntu) ]; then
                use_or_install_with_checks "VirtualBox" "VirtualBox VM manager" 'sudo apt-get -y install virtualbox virtualbox-dkms' 1
                use_or_install_with_checks "vagrant" "Vagrant VM environment manager" 'sudo apt-get -y install vagrant' 1
            else
                echo "This script is intended to be used with Ubuntu. If you are using a different Linux distro, please set up the environment manually."
                exit 1
            fi

        elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
            echo "Platform: Windows. Not currently supported."
            exit 1
        fi
    else
        echo "Aborting."
        exit 1
    fi
    install_with_blocking "vagrant plugin install vagrant-aws"
else
    echo $'Dependencies for building the viral-ngs environment are SATISFIED.\n'
    echo $'Next you have an opportunity to choose whether you would like to provision'
    echo $'the viral-ngs environment LOCALLY or on AWS EC2.'
    echo $'In either case, you will be directed to an ssh session on a VM where'
    echo $'you can work with viral-ngs.\n'

    if $(ask "Run viral-ngs LOCALLY?" Y); then
        vagrant up --provider=virtualbox
        vagrant box update
        vagrant ssh
        vagrant suspend
        if [[ $? ]]; then
            echo $'\nThe Viral-NGS virtual machine has been suspended. Next time you run this'
            echo $'script, you will be connected exactly where you left off.\n'
            exit 0
        else
            echo $'The Vial-NGS virtual machine did not close cleanly.'
            exit 1
        fi
    elif $(ask "Run viral-ngs on AWS?" Y); then
        vagrant up --provider=aws
        vagrant ssh
        exit 0
    fi
fi