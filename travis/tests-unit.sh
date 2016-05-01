#!/bin/bash
#set -e

nosetests -v \
    --logging-clear-handlers \
    --with-timer \
    --with-xunit \
    -w test/unit/
