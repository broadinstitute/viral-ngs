#!/bin/sh -f

# Derived from: https://github.com/mernst/plume-lib/blob/master/bin/trigger-travis.sh
# Which is originally under an MIT license

# Trigger a new Travis-CI job.
# Ordinarily, a new Travis job is triggered when a commit is pushed to a
# GitHub repository.  The trigger-travis.sh script provides a programmatic
# way to trigger a new Travis job.

# Usage:
#   trigger-travis.sh [--pro] [--branch BRANCH] GITHUBID GITHUBPROJECT TRAVIS_ACCESS_TOKEN [MESSAGE]
# For example:
#   trigger-travis.sh typetools checker-framework `cat ~/private/.travis-access-token` "Trigger for testing"
#
# where --pro means to use travis-ci.com instead of travis-ci.org, and
# where TRAVIS_ACCESS_TOKEN is, or ~/private/.travis-access-token contains,
# the Travis access token.
#
# Your Travis access token is the text after "Your access token is " in
# the ouput of these commands:
#   travis login && travis token
# (If the travis program isn't installed, do so with one of these two commands:
#    gem install travis
#    sudo apt-get install ruby-dev && sudo gem install travis
# Don't do "sudo apt-get install travis" which installs a trajectory analyzer.)
# Note that the Travis access token output by `travis token` differs from the
# Travis token available at https://travis-ci.org/profile .
# If you store it in in a file, make sure the file is not readable by others,
# for example by running:  chmod og-rwx ~/private

# To use this script to trigger a dependent build in Travis, do two things:
#
# 1. Set an environment variable TRAVIS_ACCESS_TOKEN by navigating to
#   https://travis-ci.org/MYGITHUBID/MYGITHUBPROJECT/settings
# The TRAVIS_ACCESS_TOKEN environment variable will be set when Travis runs
# the job, but won't be visible to anyone browsing https://travis-ci.org/.
#
# 2. Add the following before_install and after_script block to your
# .travis.yml file, where you replace OTHERGITHUB* by a specific downstream
# project, but you leave $TRAVIS_ACCESS_TOKEN as literal text:
#
# before_install:
#   - npm install --save-dev travis-after-all
#
# after_script:
#   - |
#       declare exitCode;
#       $(npm bin)/travis-after-all
#       exitCode=$?
#
#       if [ "$exitCode" -eq 0 ]; then
#         if [[ ($TRAVIS_BRANCH == master) &&
#               ($TRAVIS_PULL_REQUEST == false) ]] ; then
#           curl -LO https://raw.github.com/mernst/plume-lib/master/bin/trigger-travis.sh
#           sh trigger-travis.sh OTHERGITHUBID OTHERGITHUBPROJECT $TRAVIS_ACCESS_TOKEN
#         fi
#       fi
#
# Your .travis.yml file must not use `language: generic` because then
# npm won't be installed.
#
# Note that Travis does not fail a job if an after_success command fails.
# If you misspell a GitHub ID or project name, then this script will fail,
# but Travis won't inform you of the mistake.  So, check the end of the
# Travis build log the first time that a build succeeds.

# Here is an explanation of the conditional in the after_success block:
#
# 1. Downstream projects are triggered only for builds of the mainline, not
# branches or pull requests.  The reason is that typically a downstream
# project clones and uses the mainline.  You could enhance this script to
# accept pass an environment variable for the upstream project; the
# downstream project's build script would need to read and use that
# environment variable.  If you make this enhancement, feel free to submit
# a pull request so that others can benefit from it.
#
# 2. Downstream projects are triggered only if the Travis job number
# contains no "." or ends with ".1".  In other words, if your .travis.yml
# defines a build matrix
# (https://docs.travis-ci.com/user/customizing-the-build/#Build-Matrix)
# that runs the same job using different configurations, then the
# "after_success:" block is run only for the first configuration.
# By default an after_success: block is run for every build in the matrix,
# but you really want it to run once if all the builds in the matrix
# succeed.  Running if the first job succeeds is simple and it is usually
# adequate, even though the downstream job is triggered even if some job
# other than the first one fails.

# TODO: enable the script to clone a particular branch rather than master.
# This would require a way to know the relationships among branches in
# different GitHub projects.  It's easier to run all your tests within a
# single Travis job, if they fit within Travis's 50-minute time limit.

# An alternative to this script would be to install the Travis command-line
# client and then run:
#   travis restart -r OTHERGITHUBID/OTHERGITHUBPROJECT
# That is undesirable because it restarts an old job, destroying its history,
# rather than starting a new job which is our goal.

# Parts of this script were originally taken from
# http://docs.travis-ci.com/user/triggering-builds/

# TODO: take a --branch command-line argument.

if [ "$#" -lt 3 ] || [ "$#" -ge 7 ]; then
  echo "Wrong number of arguments $# to trigger-travis.sh; run like:"
  echo " trigger-travis.sh [--pro] [--script 'script value...'] [--branch BRANCH] GITHUBID GITHUBPROJECT TRAVIS_ACCESS_TOKEN [MESSAGE]" >&2
  exit 1
fi

if [ "$1" = "--pro" ] ; then
  TRAVIS_URL=travis-ci.com
  shift
else
  TRAVIS_URL=travis-ci.org
fi

if [ "$1" = "--branch" ] ; then
  shift
  BRANCH="$1"
  shift
else
  BRANCH=master
fi

if [ "$1" = "--script" ] ; then
  shift
  SCRIPT_VAL="$1"
  shift
else
  SCRIPT_VAL=""
fi

USER=$1
REPO=$2
TOKEN=$3
if [ $# -eq 4 ] ; then
    MESSAGE=",\"message\": \"$4\""
elif [ -n "$TRAVIS_REPO_SLUG" ] ; then
    MESSAGE=",\"message\": \"Triggered by upstream build of $TRAVIS_REPO_SLUG commit "`git rev-parse --short HEAD`"\""
else
    MESSAGE=""
fi

if [ -n "$SCRIPT_VAL" ]; then
  SCRIPT_BODY="\"script\": \"$SCRIPT_VAL\""
else
  SCRIPT_BODY=""
fi

body="{
\"request\": {
  \"branch\":\"$BRANCH\"$MESSAGE,

  \"config\":{
  $SCRIPT_BODY
}}}"

echo "Travis API request body:\n$body"

#exit 0

echo "Making request to start tests in other repository..."

# It does not work to put / in place of %2F in the URL below. 
curl -S -X POST \
  -H "Content-Type: application/json" \
  -H "Accept: application/json" \
  -H "Travis-API-Version: 3" \
  -H "Authorization: token ${TOKEN}" \
  -d "$body" \
  https://api.${TRAVIS_URL}/repo/${USER}%2F${REPO}/requests \
 | tee /tmp/travis-request-output.$$.txt

if grep -q '"@type": "error"' /tmp/travis-request-output.$$.txt; then
    exit 1
fi
if grep -q 'access denied' /tmp/travis-request-output.$$.txt; then
    exit 1
fi