nextflow.enable.dsl=2

params.rootDir = java.nio.file.Paths.get("$projectDir/../../../").toAbsolutePath().normalize().toString()
targetDir = params.rootDir + "/target/nextflow/workflows"


include { wrapper_qc } from targetDir + "/wrapper_qc/main.nf"

params.resources_test = "s3://openpipelines-bio/openpipeline_incubator/resources_test/"

workflow test_xenium {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_one",
        input: resources_test_file.resolve("spatial_qc_sample_data/xenium_tiny.qc.h5mu"),
        ingestion_method: "xenium",
        var_gene_names: "gene_ids",
        publish_dir: "test_out"
      ],
      [
        id: "sample_two",
        input: resources_test_file.resolve("spatial_qc_sample_data/xenium_tiny.qc.h5mu"),
        ingestion_method: "xenium",
        var_gene_names: "gene_ids",
        publish_dir: "test_out"
      ]
    ])

    | map{ state -> [state.id, state] }
    | wrapper_qc

    | view { output ->
        assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
        def id = output[0]
        def state = output [1]
        assert id == "combined": "Output ID should be `combined`"
        assert state instanceof Map : "State should be a map. Found: ${state}"
        assert state.containsKey("output_qc_report"): "Output should contain key `output_qc_report`"
        assert state.containsKey("output_processed_h5mu"): "Output should contain key `output_processed_h5mu`"
        assert state.output_qc_report.size() == 1 : "Expected exactly one output HTML file to be generated"
        assert state.output_qc_report.every { it.isFile()} : "All output HTML report file should exist"
        assert state.output_processed_h5mu.isDirectory() : "Output directory should exist"
        def files = state.output_processed_h5mu.listFiles().findAll { it.isFile() }
        assert files.size() == 2 : "Output directory should contain exactly 2 files, but found ${files.size()} files"
        "Output: $output"
    }
}
