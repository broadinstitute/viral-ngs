#!/bin/bash

pytest --cov-append test/integration

rc=$?; if [[ $rc != 0 ]]; then sleep 10; exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
