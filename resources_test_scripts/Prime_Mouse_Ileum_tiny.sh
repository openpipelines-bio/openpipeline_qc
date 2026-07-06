#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID="Prime_Mouse_Ileum_tiny"
OUT=resources_test/xenium/$ID

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

if [ ! -d "$OUT" ]; then
    tiny_dataset="https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.tar.gz"
    wget "$tiny_dataset" -O "$TMPDIR/xenium_tiny.tar.gz"

    mkdir -p "$TMPDIR/xenium_tiny"
    tar -xzf "$TMPDIR/xenium_tiny.tar.gz" -C "$TMPDIR/xenium_tiny"
    mkdir -p "$OUT/xenium_outs/"
    mv "$TMPDIR/xenium_tiny/Xenium_Prime_Mouse_Ileum_tiny_outs/"* "$OUT/xenium_outs/"
fi

# Create H5MU
nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial \
  -revision v0.5.0 \
  -main-script target/nextflow/convert/from_xenium_to_spatialdata//main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --input "$OUT/xenium_outs" \
  --output "$ID.zarr" \
  --publish_dir "$OUT"

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial \
  -revision v0.5.0 \
  -main-script target/nextflow/convert/from_spatialdata_to_h5mu//main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --input "$OUT/$ID.zarr" \
  --output "$ID.h5mu" \
  --publish_dir "$OUT"

# run qc workflow
nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --id "$ID" \
  --input "$OUT/$ID.h5mu" \
  --var_name_mitochondrial_genes mitochondrial \
  --var_name_ribosomal_genes ribosomal \
  --output "$ID.h5mu" \
  --publish_dir "$OUT"

find "${OUT}" -mindepth 1 ! -name "$ID.h5mu" -delete

# generate json for testing
nextflow run https://packages.viash-hub.com/vsh/openpipeline_qc \
  -r v0.3.0 \
  -main-script target/_private/nextflow/ingestion_qc/h5mu_to_qc_json/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --input "$OUT"/$ID.h5mu \
  --input "$OUT"/$ID.h5mu \
  --ingestion_method xenium \
  --min_num_nonzero_vars 1 \
  --output xenium_dataset.json \
  --output_reporting_json xenium_report_structure.json \
  --publish_dir "$OUT"

rm -f "${OUT}"/*.state.yaml

# Sync to S3
aws s3 sync \
    "$OUT" \
    s3://openpipelines-bio/openpipeline_qc/resources_test/xenium/"$ID" \
    --delete \
    --dryrun
