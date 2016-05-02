#!/bin/bash

# Modify sitecustomize.py file for coverage. Allows to cover files run in a subprocess.
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    touch "/home/travis/virtualenv/python${TRAVIS_PYTHON_VERSION}/lib/python${TRAVIS_PYTHON_VERSION}/sitecustomize.py";
fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]] && [[ $TRAVIS_PYTHON_VERSION == 2 ]]; then
    touch "/usr/local/lib/python2.7/site-packages/sitecustomize.py";
fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]] && [[ $TRAVIS_PYTHON_VERSION == 3 ]]; then
    touch "/usr/local/lib/python3.5/site-packages/sitecustomize.py";
fi

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    printf "import coverage\ncoverage.process_startup()\n" > "/home/travis/virtualenv/python${TRAVIS_PYTHON_VERSION}/lib/python${TRAVIS_PYTHON_VERSION}/sitecustomize.py";
fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]] && [[ $TRAVIS_PYTHON_VERSION == 2 ]]; then
    printf "import coverage\ncoverage.process_startup()\n" > "/usr/local/lib/python2.7/site-packages/sitecustomize.py";
fi
if [[ "$TRAVIS_OS_NAME" == "osx" ]] && [[ $TRAVIS_PYTHON_VERSION == 3 ]]; then
    printf "import coverage\ncoverage.process_startup()\n" > "/usr/local/lib/python3.5/site-packages/sitecustomize.py";
fi
