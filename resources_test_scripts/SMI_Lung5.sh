#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=SMI_Lung5
OUT=resources_test/cosmx/$ID

mkdir -p "$OUT"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/cosmx-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

if [ ! -d "$OUT" ]; then

    flat_dataset_rep_1="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep1/Lung5_Rep1+SMI+Flat+data.tar.gz"
    wget  "$flat_dataset_rep_1" -O "$TMPDIR/Lung5_Rep1.tar.gz"
    mkdir -p "$TMPDIR/Lung5_Rep1"
    tar -xzf "$TMPDIR/Lung5_Rep1.tar.gz" -C "$TMPDIR/Lung5_Rep1"
    mkdir -p "$OUT/Lung5_Rep1/"
    mv "$TMPDIR/Lung5_Rep1/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/"* "$OUT/Lung5_Rep1/"

    flat_dataset_rep_2="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    wget  "$flat_dataset_rep_2" -O "$TMPDIR/Lung5_Rep2.tar.gz"
    mkdir -p "$TMPDIR/Lung5_Rep2"
    tar -xzf "$TMPDIR/Lung5_Rep2.tar.gz" -C "$TMPDIR/Lung5_Rep2"
    mkdir -p "$OUT/Lung5_Rep2/"
    mv "$TMPDIR/Lung5_Rep2/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"* "$OUT/Lung5_Rep2/"
fi

echo "> Downloading of datasets complete"

# Subset dataset to make it tiny
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$OUT/Lung5_Rep1/"
  output : "Lung5_Rep1_tiny"
- id: Lung5_Rep2
  input: "$OUT/Lung5_Rep2/"
  output : "Lung5_Rep2_tiny"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial \
  -revision v0.5.0 \
  -main-script target/_private/nextflow/filter/subset_cosmx/main.nf \
  -params-file /tmp/params.yaml \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --publish_dir "$OUT" \
  --num_fovs 3 \
  --subset_transcripts_file True \
  --subset_polygons_file False  


echo "> Subsetting complete"

# Convert to h5mu
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$OUT/Lung5_Rep1_tiny/"
  output : "Lung5_Rep1_tiny.h5mu"
- id: Lung5_Rep2
  input: "$OUT/Lung5_Rep2_tiny/"
  output : "Lung5_Rep2_tiny.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial \
  -revision v0.5.0 \
  -main-script target/nextflow/convert/from_cosmx_to_h5mu/main.nf \
  -params-file /tmp/params.yaml \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --publish_dir "$OUT" \
  --output_compression "gzip"

echo "> Conversion to H5MU complete"

# run qc workflow
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input : "$OUT/Lung5_Rep1_tiny.h5mu"
  output: Lung5_Rep1_tiny.h5mu
- id: Lung5_Rep2
  input : "$OUT/Lung5_Rep2_tiny.h5mu"
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
publish_dir: "$OUT/"
output: Lung5_Rep2_tiny.h5mu
HERE

nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c ./src/configs/labels_ci.config


find "${OUT}" -mindepth 1 ! -name "Lung5_Rep*_tiny.h5mu" -delete

# generate json for testing
nextflow run https://packages.viash-hub.com/vsh/openpipeline_qc \
  -r v0.3.0 \
  -main-script target/_private/nextflow/ingestion_qc/h5mu_to_qc_json/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --input "$OUT"/Lung5_Rep1_tiny.h5mu \
  --input "$OUT"/Lung5_Rep2_tiny.h5mu \
  --ingestion_method cosmx \
  --min_num_nonzero_vars 1 \
  --output cosmx_dataset.json \
  --output_reporting_json cosmx_report_structure.json \
  --publish_dir "$OUT"

rm -f "${OUT}"/*.state.yaml


# Sync to S3
aws s3 sync \
    "$OUT" \
    s3://openpipelines-bio/openpipeline_qc/resources_test/cosmx/"$ID" \
    --delete \
    --dryrun
