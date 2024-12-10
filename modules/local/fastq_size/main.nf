process FASTQ_SIZE {
    
    tag "$meta.id"
    label 'process_low'

    // Define the input channel
    input:
    tuple val(meta), path(file_in)

    // Define the output channel
    output:
    tuple val(meta), path(file_in), stdout , emit: outs

    // Define the command to run
    script:
    """
    gzip -l ${file_in} | awk 'NR==2 {print \$2}'
    """
}