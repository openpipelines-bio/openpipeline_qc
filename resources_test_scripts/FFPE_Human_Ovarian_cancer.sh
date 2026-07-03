#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID="FFPE_Human_Ovarian_Cancer_tiny"
OUT=resources_test/visium/$ID


# Input Files - download to the specific directory
curl -o "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -o "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -o "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Extract in the specific directory
tar xvf "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" -C "$OUT/visium_outs/"

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir -p "$OUT/visium_tiny/Visium_FFPE_Human_Ovarian_Cancer_tiny"
convert "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" -resize 2000x2000 "$OUT/visium_tiny/Visium_FFPE_Human_Ovarian_Cancer_image_tiny.jpg"
for f in "$OUT"/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do 
  gzip -cdf "$f" | head -n 40000 | gzip -c > "$OUT/visium_tiny/Visium_FFPE_Human_Ovarian_Cancer_tiny/$(basename "$f")"; 
done

echo "> Downloading and subsampling of datasets complete"

# Run spaceranger
nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial \
  -revision v0.3.0 \
  -profile docker \
  -resume \
  -c src/configs/labels_ci.config \
  -main-script target/nextflow/workflows/ingestion/spaceranger_mapping/main.nf \
  --id $ID \
  --input "$OUT/visium_tiny/Visium_FFPE_Human_Ovarian_Cancer_tiny" \
  --gex_reference "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" \
  --probe_set "$OUT/visium_outs/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" \
  --image "$OUT/visium_tiny/Visium_FFPE_Human_Ovarian_Cancer_image_tiny.jpg" \
  --slide "V10L13-020" \
  --area "D1" \
  --create_bam "false" \
  --output_raw "$ID" \
  --output_h5mu "$ID.h5mu" \
  --publish_dir "$OUT"

# run qc workflow
cat > /tmp/params.yaml << HERE
id: $ID
input : "$OUT/$ID.h5mu"
output: $ID.h5mu
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
publish_dir: "$OUT"
HERE

nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c ./src/configs/labels_ci.config

find "${OUT}" -mindepth 1 ! -name "$ID.h5mu" -delete

# generate json for testing
nextflow run https://packages.viash-hub.com/vsh/openpipeline_qc \
  -latest \
  -r v0.3.0 \
  -main-script target/_private/nextflow/ingestion_qc/h5mu_to_qc_json/main.nf \
  -c src/configs/labels_ci.config \
  -profile docker \
  -resume \
  --input "$OUT"/$ID.h5mu \
  --input "$OUT"/$ID.h5mu \
  --ingestion_method visium \
  --min_num_nonzero_vars 1 \
  --output visium_dataset.json \
  --output_reporting_json visium_report_structure.json \
  --publish_dir "$OUT"

rm -f "${OUT}"/*.state.yaml

# Sync to S3
aws s3 sync \
    --profile di \
    "$OUT" \
    s3://openpipelines-bio/openpipeline_qc/resources_test/visium/"$ID" \
    --delete \
    --dryrun
