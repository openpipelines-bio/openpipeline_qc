ABSOLUTE_INPUT_DATA=$(realpath $par_input_data)
ABSOLUTE_INPUT_STRUCTURE=$(realpath $par_input_structure)
ABSOLUTE_OUTPUT=$(realpath $par_output_qc_report)

# The report is built inside the app's folder, so render from a copy of its
# own: the one in the image belongs to root, and is read-only under Apptainer.
SIQC_DIR=$(mktemp -d "$meta_temp_dir/siqc_XXXXXX")
trap 'rm -rf "$SIQC_DIR"' EXIT
cp -r /opt/siqc/. "$SIQC_DIR"
cd "$SIQC_DIR"
mkdir -p src/data

npm run cli render -- --data "$ABSOLUTE_INPUT_DATA" --structure "$ABSOLUTE_INPUT_STRUCTURE" --output "$ABSOLUTE_OUTPUT"
