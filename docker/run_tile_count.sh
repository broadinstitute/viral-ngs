#!/bin/bash
set -e -o pipefail

if [ $# -eq 0 ]; then
    echo "Usage: $(basename $0) path/to/RunInfo.xml"
    exit 0
fi

if [ ! -e "$1" ]; then
    echo "The specified file does not exist: $1"
    exit 1
fi

# Parse the lane count & run ID from RunInfo.xml file
lane_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@LaneCount)" $1)
if [ -z "$lane_count" ]; then
    echo "Could not parse LaneCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
fi

surface_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SurfaceCount)" $1)
if [ -z "$surface_count" ]; then
    echo "Could not parse SurfaceCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
fi

swath_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SwathCount)" $1)
if [ -z "$swath_count" ]; then
    echo "Could not parse SwathCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
fi

tile_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@TileCount)" $1)
if [ -z "$tile_count" ]; then
    echo "Could not parse TileCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
fi

# total data size more roughly tracks total tile count
total_tile_count=$((lane_count*surface_count*swath_count*tile_count))

echo $total_tile_count