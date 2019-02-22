#!/bin/bash

set -eu -o pipefail

pytest --cov-append $PYTEST_EXTRA_OPTS test

rc=$?; if [[ $rc != 0 ]]; then sleep 10; exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
