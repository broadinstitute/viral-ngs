#!/bin/bash

#
# Calculate the memory allocated to the process in mb and gb.
#
# Usage: calc_mem.sh MEM_UNIT MEM_FRACTION
#
# MEM_UNIT is mb or gb
# MEM_FRACTION is an integer 1-100
#

set -eu -o pipefail 

if [ "$#" -ne 2 ]; then
    echo "Usage: calc_mem.sh MEM_UNIT MEM_FRACTION"
    exit 1
fi

MEM_UNIT=$1
MEM_FRACTION=$2

if [ "$MEM_UNIT" = "mb" ]
then
    DIVIDER=1
elif [ "$MEM_UNIT" = "gb" ]
then
    DIVIDER=1024
else
    echo "Invalid unit: $MEM_UNIT"
    exit 1
fi

if [ -e /sys/fs/cgroup/memory/memory.limit_in_bytes ]
then
    MEM_IN_BYTES=`cat /sys/fs/cgroup/memory/memory.limit_in_bytes`
else
    MEM_IN_BYTES=$(head -n1 /proc/meminfo | awk '{printf "%d", int($2*1024)}' )
fi

echo $MEM_IN_BYTES $MEM_FRACTION $DIVIDER | awk '{printf "%d", int($1*$2/100/1024/1024/$3)}'
