#!/bin/bash

set -eo pipefail


# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_anticmv
OUT=resources_test/cellranger/$ID

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# dataset page:
# https://www.10xgenomics.com/resources/datasets/integrated-gex-totalseqc-and-tcr-analysis-of-connect-generated-library-from-5k-cmv-t-cells-2-standard

# check whether reference is available
reference_dir="resources_test/cellranger/ref_gencodev41_chr1"
genome_tar="$reference_dir/reference_cellranger.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
    exit 1
fi

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/5k_human_antiCMV_T_TBNK_connect_Multiplex"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_fastqs.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit head -n 200000 "$input" | gzip > "$output"
  fi
}

orig_sample_id="5k_human_antiCMV_T_TBNK_connect"

seqkit_head "$tar_dir/gex_1/${orig_sample_id}_GEX_1_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_GEX_1_subset_S1_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/gex_1/${orig_sample_id}_GEX_1_S1_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_GEX_1_subset_S1_L001_R2_001.fastq.gz"

# download immune panel fasta if needed
feature_reference="$raw_dir/feature_reference.csv"
if [[ ! -f "$feature_reference" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv" -O "$feature_reference"
fi

# Run mapping pipeline
nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/workflows/ingestion/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --id "$ID" \
  --input "$raw_dir" \
  --library_id "${orig_sample_id}_GEX_1_subset" \
  --library_type "Gene Expression" \
  --gex_reference "$genome_tar" \
  --feature_reference "$feature_reference" \
  --output_h5mu "${orig_sample_id}.h5mu" \
  --publish_dir "$OUT"

rm -rf "${OUT}/${ID}.cellranger_multi.output_raw/"


# run qc workflow
nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --id "$ID" \
  --input "$OUT/$orig_sample_id.h5mu" \
  --var_name_mitochondrial_genes mitochondrial \
  --var_name_ribosomal_genes ribosomal \
  --output "${orig_sample_id}_qc.h5mu" \
  --publish_dir "$OUT"

# run cellbender
nextflow \
  run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/correction/cellbender_remove_background/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --id "$ID" \
  --input "${OUT}/${orig_sample_id}_qc.h5mu" \
  --output "${orig_sample_id}_qc_cellbender.h5mu" \
  --epochs 5 \
  --output_compression gzip \
  --publish_dir "$OUT"

# Subset h5mu
nextflow run https://packages.viash-hub.com/vsh/openpipeline \
  -r v4.1.1 \
  -main-script target/nextflow/filter/subset_h5mu/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --id "${orig_sample_id}_10k" \
  --input "${OUT}/${orig_sample_id}_qc_cellbender.h5mu" \
  --number_of_observations 10000 \
  --output "${orig_sample_id}_10k.h5mu" \
  --output_compression gzip \
  --publish_dir "$OUT"

find "${OUT}" -mindepth 1 ! -name "${orig_sample_id}_10k.h5mu" -delete

cat > /tmp/add_metadata_obs.py <<EOF
import mudata as mu
import numpy as np
import pandas as pd
import os

# List of h5mu files
h5mu_file = "${OUT}/${orig_sample_id}_10k.h5mu"

# Metadata values to randomly assign
donor_ids = ["donor_1", "donor_2", "donor_3"]
cell_types = ["CD4+ T cell", "CD8+ T cell", "B cell", "NK cell", "Monocyte"]
batches = ["batch_A", "batch_B"]
conditions = ["treated", "control"]
    
# Load MuData object
mdata = mu.read_h5mu(h5mu_file)
rna = mdata.mod["rna"]
n_obs = rna.n_obs

# Generate random metadata
np.random.seed(42)

# Create metadata
rna.obs["donor_id"] = np.random.choice(donor_ids, size=n_obs)
rna.obs["cell_type"] = np.random.choice(cell_types, size=n_obs)
rna.obs["batch"] = np.random.choice(batches, size=n_obs)
rna.obs["condition"] = np.random.choice(conditions, size=n_obs)

# Add a continuous variable too
rna.obs["quality_score"] = np.random.uniform(0, 1, size=n_obs)

# Save the modified MuData object
mu.write_h5mu(h5mu_file, mdata)
print(f"Added metadata to {h5mu_file}")

print("All files processed successfully!")
EOF

python /tmp/add_metadata_obs.py

# generate json for testing
nextflow run https://packages.viash-hub.com/vsh/openpipeline_qc \
  -r v0.3.0 \
  -main-script target/_private/nextflow/ingestion_qc/h5mu_to_qc_json/main.nf \
  -resume \
  -profile docker,mount_temp \
  -c ./src/configs/labels_ci.config \
  --input "${OUT}/${orig_sample_id}_10k.h5mu" \
  --input "${OUT}/${orig_sample_id}_10k.h5mu" \
  --ingestion_method cellranger_multi \
  --obs_metadata "donor_id;cell_type;batch;condition" \
  --output sc_dataset.json \
  --output_reporting_json sc_report_structure.json \
  --publish_dir "$OUT"

rm -f "${OUT}"/*.state.yaml

aws s3 sync \
  "$OUT" \
  s3://openpipelines-bio/openpipeline_qc/resources_test/cellranger/"$ID" \
  --delete \
  --dryrun
