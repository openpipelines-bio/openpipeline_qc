# openpipeline_qc 0.2.0

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