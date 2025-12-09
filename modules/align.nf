// modules/align.nf - FASTQ alignment and variant calling using Dragen

process ALIGN {
    tag "${sample_id}"
    label 'dragen'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(panel), path(fastq_list)

    output:
    tuple val(sample_id), val(panel), path("${sample_id}.cram"), path("${sample_id}.cram.crai"), emit: cram
    path "${sample_id}.hard-filtered.gvcf.gz*", emit: gvcf
    path "${sample_id}*.csv", emit: metrics
    path "${sample_id}.hla.*", emit: hla, optional: true

    script:
    def panel_config = params.panels[panel]
    def bed_file = panel_config.bed
    def human_ref = panel_config.human_ref
    def staging_dir = params.dragen_staging_dir ?: '/staging/tmp'

    // Panel-specific options
    def umi_opts = panel == 'HEME' ? """
        --umi-enable true \\
        --umi-source qname \\
        --umi-correction-scheme random \\
        --umi-min-supporting-reads 1 \\
        --umi-metrics-interval-file ${bed_file} \\
        --vc-enable-umi-germline true \\
        --vc-enable-high-sensitivity-mode true""" : """
        --enable-duplicate-marking true"""

    """

    dragen \\
        --intermediate-results-dir ${staging_dir} \\
        --enable-map-align true \\
        --enable-map-align-output true \\
        --output-format CRAM \\
        --generate-sa-tags true \\
        --enable-sort true \\
        --soft-read-trimmers polyg,quality \\
        --trim-min-quality 2 \\
        --ref-dir ${human_ref} \\
        --qc-coverage-tag-1 target_bed \\
        --qc-coverage-region-1 ${bed_file} \\
        --qc-coverage-reports-1 cov_report \\
        --qc-coverage-ignore-overlaps true \\
        --enable-variant-caller true \\
        --vc-combine-phased-variants-distance 6 \\
        --vc-emit-ref-confidence GVCF \\
        --enable-hla true \\
        --fastq-list ${fastq_list} \\
        --fastq-list-sample-id ${sample_id} \\
        --output-directory . \\
        --output-file-prefix ${sample_id} \\
        ${umi_opts}

    echo "[\$(date)] Alignment complete for ${sample_id}"
    """
}
