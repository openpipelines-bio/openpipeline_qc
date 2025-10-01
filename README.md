# OpenPipeline QC

OpenPipeline QC contains an extensible workflow for the generation of a stand-alone, interactive ingestion QC report for single-cell and spatial omics data. It is part of the [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) ecosystem and extends it with specialized functionality for ingestion quality control.

[![ViashHub](https://img.shields.io/badge/ViashHub-openpipeline_qc-7a4baa.svg)](https://www.viash-hub.com/packages/openpipeline_qc)
[![GitHub](https://img.shields.io/badge/GitHub-openpipelines--bio%2Fopenpipeline_qc-blue.svg)](https://github.com/openpipelines-bio/openpipeline_qc)
[![GitHub License](https://img.shields.io/github/license/openpipelines-bio/openpipeline_qc.svg)](https://github.com/openpipelines-bio/openpipeline_qc/blob/main/LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/openpipelines-bio/openpipeline_qc.svg)](https://github.com/openpipelines-bio/openpipeline_qc/issues)
[![Viash version](https://img.shields.io/badge/Viash-v0.9.3-blue.svg)](https://viash.io)

## Functionality

### 1. QC Metrics Calculation

The workflow calculates extensive ingestion QC metrics, including:

**Observation level metrics:**
- `total_{var_type}_by_{expr_type}` - E.g. "total_genes_by_counts". Number of genes with positive counts in a cell.
- `total_{expr_type}` - E.g. "total_counts". Total number of counts for a cell.
- `pct_{expr_type}_in_top_{n}_{var_type}` - E.g. "pct_counts_in_top_50_genes". Cumulative percentage of counts for 50 most expressed genes in a cell.
- `total_{expr_type}_{qc_var}` - E.g. "total_counts_mito". Total number of counts for variables in qc_vars.
- `pct_{expr_type}_{qc_var}` - E.g. "pct_counts_mito". Proportion of total counts for a cell which are mitochondrial.

**Variable level metrics:**
- `total_{expr_type}` - E.g. "total_counts". Sum of counts for a gene.
- `n_genes_by_{expr_type}` - E.g. "n_genes_by_counts". The number of genes with at least 1 count in a cell.
- `mean_{expr_type}` - E.g. "mean_counts". Mean expression over all cells.
- `n_cells_by_{expr_type}` - E.g. "n_cells_by_counts". Number of cells this expression is measured in.
- `pct_dropout_by_{expr_type}` - E.g. "pct_dropout_by_counts". Percentage of cells this feature does not appear in.

**Spatial data metrics:**
For spatial data, the workflow also surfaces segmentation metrics.

### 2. Interactive Report Generation

After calculating QC metrics, the workflow generates a standalone interactive report using the [siqc](https://github.com/openpipelines-bio/siqc) package (standalone, interactive QC reporting engine). 

**siqc features:**
- 🚀 **Progressive Loading**: Large datasets (millions of cells) load incrementally without blocking the UI
- 📊 **Interactive Visualizations**: Histograms, scatter plots, heatmaps, and spatial plots with Plotly.js
- 🗜️ **Efficient Compression**: Binary data format with gzip compression for optimal file sizes
- 📄 **Standalone Reports**: Single HTML file with no external dependencies
- 🎯 **Spatial Analysis**: Advanced spatial plotting with faceting and optimization
- ⚡ **Performance Optimized**: Web Workers, typed arrays, and multi-scale rendering

### 3. Quality Assessment for Downstream Processing

The generated report can be used to set thresholds for QC-based filtering, which can then be applied for further processing of the data using the core [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) of [OpenPipeline Spatial](https://github.com/openpipelines-bio/openpipeline_spatial/) packages.

**Supported platforms:**
- 10x Genomics single-cell multi-omics (Cellranger Multi output)
- 10x Genomics Xenium spatial transcriptomics

## Execution via CLI or Seqera Cloud

The openpipeline_qc package is available via [Viash Hub](https://www.viash-hub.com/packages/openpipeline_qc/latest/), where you can receive instructions on how to run the end-to-end workflow as well as individual subworkflows or components.

It's possible to run the workflow directly from Seqera Cloud. The necessary Nextflow schema files have been [built and provided with the workflows](https://packages.viash-hub.com/vsh/openpipeline_qc/-/tree/build/main/target/nextflow/workflows/generate_qc_report) in order to use the form-based input. However, Seqera Cloud can not deal with multiple-value parameters for batch processing of multiple samples. Therefore, it's better to use Viash Hub also here for launching the workflow on Seqera Cloud.

* Navigate to the [Viash Hub package page](https://www.viash-hub.com/packages/openpipeline_qc/latest/), select the workflow you want to launch and click the `launch` button.
* Select the execution environment of choice (e.g. `Seqera Cloud`, `Nextflow` or `Executable`)
* Fill in the form with the required parameters and launch the workflow.

## Support

For issues specific to quality control analysis, please use the [GitHub issues tracker](https://github.com/openpipelines-bio/openpipeline_qc/issues). For general OpenPipeline questions, refer to the main [OpenPipeline documentation](https://openpipelines.bio/).
