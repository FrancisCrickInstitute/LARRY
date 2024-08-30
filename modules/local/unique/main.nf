process UNIQUE {
    
    tag "$meta.id"
    label 'process_low'

    // Define the input channel
    input:
    tuple val(meta), path(file_in)

    // Define the output channel
    output:
    tuple val(meta), path("*unique*"), emit: file_out

    // Define the command to run
    script:
    def args        = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_name        = file_in.toString()
    def out_name         = "${prefix}_unique_${file_name}"

    """
    sort ${file_name} | uniq | sed 's/-1//g' > ${out_name}
    """
}