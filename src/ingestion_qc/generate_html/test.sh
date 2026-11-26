echo ">> Generating report"
"$meta_executable" \
    --input_data "$meta_resources_dir/sc_dataset.json" \
    --input_structure "$meta_resources_dir/sc_report_structure.json" \
    --output_qc_report "index.html" \

echo ">> Checking output"
[ ! -f "index.html" ] && echo "Error: Output report does not exist." && exit 1

# The report app is installed as root when the image is built. A container
# that runs as another user, as Docker's `-u` and rootless Podman do, must
# still be able to render; `node` is the image's own non-root user.
echo ">> Generating report as a user other than root"
chmod a+rwx .
su node --shell /bin/bash --command "\"$meta_executable\" \
    --input_data \"$meta_resources_dir/sc_dataset.json\" \
    --input_structure \"$meta_resources_dir/sc_report_structure.json\" \
    --output_qc_report \"$PWD/index_not_root.html\""

echo ">> Checking output"
[ ! -f "index_not_root.html" ] && echo "Error: Output report does not exist when rendered as another user." && exit 1

echo ">> Test succesful" && exit 0
