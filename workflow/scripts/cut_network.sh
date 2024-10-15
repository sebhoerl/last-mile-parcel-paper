#!/bin/bash
set -e

# This script uses osmosis to cut the OpenSteetMap data

# Show version to test if osmosis is executable
osmosis -v

# Extract road network
osmosis --read-pbf ${snakemake_input[osm]} --tag-filter accept-ways highway=* --bounding-polygon file=${snakemake_input[polygon]} completeWays=yes --used-node --write-xml ${snakemake_output[0]}
