#!/usr/bin/env bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# build the component
viash ns build --setup cb -q wrapper_qc

# run the test
nextflow run . \
  -main-script src/workflows/wrapper_qc/test.nf \
  -profile docker,no_publish,local \
  -entry test_xenium \
  -c src/configs/labels_ci.config \
  -resume
