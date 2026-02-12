workflow run_wf {
  take: input_ch
  main:
    output_ch = input_ch
    | generate_qc_report.run(
        fromState: [
            id: "id",
            input: "input",
            ingestion_method: "ingestion_method",
            sample_metadata: "sample_metadata",
            var_gene_names: "var_gene_names",
            obs_metadata: "obs_metadata"
        ],
        toState: [
            output_qc_report: "output_qc_report",
            output_processed_h5mu: "output_processed_h5mu"
        ]
    )
  emit: output_ch
}
