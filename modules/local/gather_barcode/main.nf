process GATHER_BARCODE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq) , val(ham_larry) , val(within_cell_cutoff), val(within_clone_cutoff), val(min_larry_umi)

    output:
    tuple val(meta), path("*clone_output.csv"), emit: outs

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template "LARRY_barcode_preprocessing_nextflow.py"

}