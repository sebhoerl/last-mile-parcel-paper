#!/bin/bash
set -e 

## This script builds VROOM with our custom additions

# remove directory created by snakemake
rm -rf results/vroom

# clone
git clone https://github.com/sebhoerl/vroom.git results/vroom

# checkout specific branch
cd results/vroom
git checkout c65d0a39

# initialize submodules
cd src
git submodule init
git submodule update

# build
make
