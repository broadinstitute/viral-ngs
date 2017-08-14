#!/bin/bash
set -e

echo "pip installing required python packages"
#pip install -r requirements.txt

#PYVER=`python -V 2>&1 | cut -c 8`
PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`
if [ "$PYVER" = "3" ]; then
    echo "pip installing snakemake packages (py3 only)"
    pip install --quiet -r requirements-py3.txt
elif [ "$PYVER" = "2" ]; then
    echo "pip install py2 packages"
    pip install --quiet -r requirements-py2.txt
fi

python --version

echo "pip installing test-related packages (coveralls, etc.)"
pip install --quiet -r requirements-py-tests.txt
