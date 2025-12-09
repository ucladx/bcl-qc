// modules/demux.nf - BCL to FASTQ demultiplexing using Dragen

process DEMUX {
    tag "${index}"
    label 'dragen'

    publishDir "${params.fastqs_dir}/${run_dir.name}/", mode: 'copy'

    input:
    path run_dir
    tuple val(index), val(panel), path(samplesheet)

    output:
    tuple val(index), val(panel), path("fastq_list.csv"), emit: fastq_list
    path "Reports/*", emit: reports
    path "*.fastq.gz", emit: fastqs

    when:
    panel != 'UNKNOWN'

    script:
    """

    dragen \\
        --bcl-conversion-only true \\
        --bcl-use-hw false \\
        --bcl-only-matched-reads true \\
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --output-directory ${index}

    echo "[\$(date)] Demux complete"
    """
}
