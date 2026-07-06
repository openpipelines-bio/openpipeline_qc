# openpipeline_qc (unreleased)

## MINOR CHANGES

* Migration of test resources to the package-specific `s3://openpipelines-bio/openpipeline_qc/resources_test` bucket (PR #26):

  - Add `.info.test_resources` to `_viash.yaml` to specify where test resources need to be synced from.
  - Test resources were regenerated via the scripts in `resources_test_scripts/`.

# openpipeline_qc 0.3.0

## MAJOR CHANGES

* Bump OpenPipeline dependency to v4.1.0, relevant updates include major changes to memory consumption and runtimes for and support for MuData encoded in Zarr format for `calculate_qc_metrics` (PR #24). 

* Bump CraftBox dependency to v0.3.2, relevant updates include support for symbolic link preservation for `move_files_to_directory` (PR #24).

* Bump Viash to 0.9.7 (PR #24).

* Bump anndata to 0.12.16 (PR #24).

* Bump mudata to 0.3.8 (PR #24).

* Testing: bump viashpy to 0.10.0 (PR #24).

# openpipeline_qc 0.2.2

## BUG FIX

* `src/ingestion_qc/h5mu_to_qc_json`: Use float32 instead of float16 for CosMx staining metrics to prevent overflow to infinity (PR #17).

# openpipeline_qc 0.2.0

## NEW FEATURES

* Add support for CosMx data in the QC Reporting workflow (PR #12).

## MINOR CHANGES

* Added support for Visium data ingestion and QC reporting (PR #16).

# openpipeline_qc 0.2.0

## MINOR CHANGES

* `src/workflows/generate_qc_report`: Add min and max constraints to `--max_samples_per_report`, to limit the number of samples per report between 1 and 20 (PR #11).

* Update `openpipeline` dependency to v4.0.0 and `siqc` to version v0.2.0 (PR #15).

# openpipeline_qc 0.1.1

## MINOR CHANGES

* Add a README (PR #9).

# openpipeline_qc 0.1.0

Initial release containing a QC Reporting workflow for Xenium or CellRanger Multi ingested data, with corresponding components.