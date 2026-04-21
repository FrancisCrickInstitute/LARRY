process GATHER_BARCODE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), file("*.png"), path("*_clone_output.csv"), emit: outs
    tuple val(meta), path("*_qc_exclusions.csv"), path("*_qc_summary.csv"), emit: qc

    script:

    template "LARRY_barcode_preprocessing_nextflow.py"

}