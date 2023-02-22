#!/bin/bash
set -eu -o pipefail

#
# Script: tests-docker.sh
#
# Run tests on the docker image
#

# Tests of util.misc.available_cpu_count()
for TEST_CPUS in $(seq 1 2)
do		 
    REPORTED_CPUS=$(docker run --cpus $TEST_CPUS local/build-container:build /bin/bash -c \
			    'cd source && python -c "import util.misc; print(util.misc.available_cpu_count())"')
    if [[ $REPORTED_CPUS -ne $TEST_CPUS ]]
    then
        echo "Problem with util.misc.available_cpu_count: reports $REPORTED_CPUS instead of $TEST_CPUS"
	exit 1
    else
	echo "util.misc.available_cpu_count correctly reported $REPORTED_CPUS cpus."
    fi
done
