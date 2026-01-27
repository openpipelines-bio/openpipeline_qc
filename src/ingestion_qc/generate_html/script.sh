ABSOLUTE_INPUT_DATA=$(realpath $par_input_data)
ABSOLUTE_INPUT_STRUCTURE=$(realpath $par_input_structure)
ABSOLUTE_OUTPUT=$(realpath $par_output_qc_report)

cd /opt/siqc
mkdir src/data

npm run cli render -- --data "$ABSOLUTE_INPUT_DATA" --structure "$ABSOLUTE_INPUT_STRUCTURE" --output "$ABSOLUTE_OUTPUT"
