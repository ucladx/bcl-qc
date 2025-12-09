// modules/multiqc.nf - Aggregate QC reports with MultiQC

process MULTIQC {
    tag "multiqc"
    label 'low_memory'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fastqs_dir
    path bams_dir
    path qcsum_files
    path interop_plot

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    def qcsum_csv = qcsum_files ? "cat ${qcsum_files.join(' ')} | head -1 > qcsum_mqc.csv && tail -n +2 -q ${qcsum_files.join(' ')} >> qcsum_mqc.csv" : "echo 'No qcsum files'"
    """
    # Remove WGS coverage reports (irrelevant for targeted panels)
    find ${bams_dir} -name "*.wgs_*.csv" -delete 2>/dev/null || true

    # Aggregate qcsum files if available
    if [ -n "${qcsum_files}" ] && [ "${qcsum_files}" != "[]" ]; then
        # Get header from first file
        head -1 \$(echo ${qcsum_files} | tr ' ' '\\n' | head -1) > qcsum_mqc.csv
        # Append data from all files (skip headers)
        for f in ${qcsum_files}; do
            tail -n +2 "\$f" >> qcsum_mqc.csv
        done
    fi

    multiqc \\
        --force \\
        --config ${projectDir}/config/multiqc_config.yaml \\
        --outdir . \\
        ${bams_dir} \\
        ${fastqs_dir}
    """
}
