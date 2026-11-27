# The report app is installed as root when the image is built. A container
# that runs as another user, as Docker's `-u` and rootless Podman do, must
# still be able to render; `node` is the image's own non-root user. It renders
# first: a render as root would leave folders behind that only root can write
# to. The test folder is mounted from the host, where other users need not be
# able to read it, so the component and its inputs are copied to a folder
# `node` can use.
echo ">> Generating report as a user other than root"
NOT_ROOT_DIR=$(mktemp -d)
cp "$meta_executable" \
    "$meta_resources_dir/sc_dataset.json" \
    "$meta_resources_dir/sc_report_structure.json" \
    "$NOT_ROOT_DIR"
chmod -R a+rwx "$NOT_ROOT_DIR"
su node --shell /bin/bash --command "cd \"$NOT_ROOT_DIR\" && ./$(basename "$meta_executable") \
    --input_data sc_dataset.json \
    --input_structure sc_report_structure.json \
    --output_qc_report index_not_root.html"

echo ">> Checking output"
[ ! -f "$NOT_ROOT_DIR/index_not_root.html" ] && echo "Error: Output report does not exist when rendered as another user." && exit 1

echo ">> Generating report"
"$meta_executable" \
    --input_data "$meta_resources_dir/sc_dataset.json" \
    --input_structure "$meta_resources_dir/sc_report_structure.json" \
    --output_qc_report "index.html" \

echo ">> Checking output"
[ ! -f "index.html" ] && echo "Error: Output report does not exist." && exit 1

echo ">> Test succesful" && exit 0
