process GATHER_BARCODE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*clone_output.csv"), emit: outs

    script:

    template "LARRY_barcode_preprocessing_nextflow.py"

}