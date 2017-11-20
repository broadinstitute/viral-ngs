#!/bin/bash

pytest test/integration

rc=$?; if [[ $rc != 0 ]]; then sleep 4 && exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
