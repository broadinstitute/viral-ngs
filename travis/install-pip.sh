#!/bin/sh

echo "pip installing required python packages"
pip install -r requirements.txt

PYVER=`python --version | cut -c 8`
if [ $PYVER == "3" ]; then
    echo "pip installing snakemake packages (py3 only)"
    pip install -r requirements-pipes.txt
fi

echo "pip installing coveralls"
pip install -q coveralls nose-cov
