#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ============================================================================
// BCL-QC Pipeline
// ============================================================================
// Demultiplexing, alignment, and QC for Illumina NovaSeq clinical NGS data
// ============================================================================

log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║                      BCL-QC Pipeline                         ║
    ╚══════════════════════════════════════════════════════════════╝
    run_dir     : ${params.run_dir}
    sampleinfo  : ${params.sampleinfo ?: 'not provided'}
    fastqs_dir  : ${params.fastqs_dir}
    bams_dir    : ${params.bams_dir}
    ─────────────────────────────────────────────────────────────────
    """.stripIndent()

// ----------------------------------------------------------------------------
// Input validation
// ----------------------------------------------------------------------------

def validateInputs() {
    def errors = []

    // Required: run_dir
    if (!params.run_dir) {
        errors << "ERROR: --run_dir is required"
    } else {
        def run_path = file(params.run_dir)
        if (!run_path.exists()) {
            errors << "ERROR: run_dir does not exist: ${params.run_dir}"
        } else if (!run_path.isDirectory()) {
            errors << "ERROR: run_dir is not a directory: ${params.run_dir}"
        }
    }

    // Validate sampleinfo if provided
    if (params.sampleinfo) {
        def si = file(params.sampleinfo)
        if (!si.exists()) {
            errors << "ERROR: sampleinfo file does not exist: ${params.sampleinfo}"
        }
    }

    // sampleinfo required for QC step
    if (!params.skip_qc && !params.sampleinfo) {
        log.warn "WARNING: --sampleinfo not provided. QC step will skip per-sample metrics."
    }

    // Check reference files exist
    def picard_ref = file(params.picard_ref)
    if (!params.skip_qc && !picard_ref.exists()) {
        errors << "ERROR: Picard reference not found: ${params.picard_ref}"
    }

    // Validate panel configurations
    ['PCP', 'HEME'].each { panel ->
        def cfg = params.panels[panel]
        if (cfg) {
            def bed = file(cfg.bed)
            def ref = file(cfg.human_ref)
            if (!bed.exists()) {
                log.warn "WARNING: ${panel} BED file not found: ${cfg.bed}"
            }
            if (!ref.exists()) {
                log.warn "WARNING: ${panel} reference not found: ${cfg.human_ref}"
            }
        }
    }

    // Check for samplesheets in run_dir (unless skipping demux)
    if (!params.skip_demux && params.run_dir) {
        def run_path = file(params.run_dir)
        if (run_path.exists()) {
            def samplesheets = run_path.listFiles().findAll { it.name =~ /SampleSheet_.*\.csv/ }
            if (samplesheets.isEmpty()) {
                errors << "ERROR: No SampleSheet_*.csv files found in run_dir: ${params.run_dir}"
            }
        }
    }

    if (errors) {
        errors.each { log.error it }
        error "Input validation failed. See errors above."
    }
}

validateInputs()

run_dir = file(params.run_dir, checkIfExists: true)
run_name = run_dir.name

// Set output directory
def final_outdir = params.outdir ?: "${params.bams_dir}/${run_name}"

// ----------------------------------------------------------------------------
// Include modules
// ----------------------------------------------------------------------------

include { DEMUX }           from './modules/demux'
include { ALIGN }           from './modules/align'
include { PICARD_HS }       from './modules/picard_hs'
include { QCSUM_METRICS }   from './modules/qcsum_metrics'
include { INTEROP_PLOT }    from './modules/interop_plot'
include { MULTIQC }         from './modules/multiqc'

// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

def get_panel_from_index(index) {
    // Determine panel type from samplesheet index naming convention
    if (index.contains('U5N2_I10')) {
        return 'HEME'
    } else if (index.contains('I10')) {
        return 'PCP'
    } else if (index.contains('I8N2_N10')) {
        return 'EXOME'  // Will be skipped in alignment
    }
    return 'UNKNOWN'
}

def get_index_from_samplesheet(samplesheet) {
    // Extract index from SampleSheet_{index}.csv filename
    def name = samplesheet.name
    def matcher = name =~ /SampleSheet_(.+)\.csv/
    return matcher ? matcher[0][1] : null
}

// ----------------------------------------------------------------------------
// Main workflow
// ----------------------------------------------------------------------------

workflow {

    // ========================================================================
    // DEMUX: BCL to FASTQ conversion
    // ========================================================================

    if (!params.skip_demux) {
        // Find all samplesheets in run directory
        samplesheets_ch = Channel
            .fromPath("${params.run_dir}/SampleSheet_*.csv")
            .map { ss ->
                def index = get_index_from_samplesheet(ss)
                def panel = get_panel_from_index(index)
                [index, panel, ss]
            }

        DEMUX(
            run_dir,
            samplesheets_ch
        )

        fastq_lists_ch = DEMUX.out.fastq_list
    } else {
        // Use existing FASTQs
        fastq_lists_ch = Channel
            .fromPath("${params.fastqs_dir}/${run_name}/*/fastq_list.csv")
            .map { fq_list ->
                def index = fq_list.parent.name
                def panel = get_panel_from_index(index)
                [index, panel, fq_list]
            }
    }

    // ========================================================================
    // ALIGN: FASTQ to CRAM with variant calling
    // ========================================================================

    if (!params.skip_align) {
        // Parse fastq_list.csv to get individual samples
        samples_ch = fastq_lists_ch
            .filter { index, panel, fq_list -> panel != 'EXOME' }  // Skip exome
            .flatMap { index, panel, fq_list ->
                def samples = fq_list.text
                    .readLines()
                    .drop(1)  // Skip header
                    .collect { it.split(',')[1] }  // RGSM column
                    .unique()
                samples.collect { sample_id -> [sample_id, panel, fq_list] }
            }

        ALIGN(samples_ch)

        aligned_ch = ALIGN.out.cram
    } else {
        // Use existing BAMs/CRAMs from sampleinfo
        if (params.sampleinfo) {
            aligned_ch = Channel
                .fromPath(params.sampleinfo)
                .splitCsv(sep: '\t', header: true)
                .map { row ->
                    def sample_id = row['Samples']
                    def panel = row['Panel'] ?: row['Tumor']
                    def bam_path = row['BAM Path']
                    [sample_id, panel, file(bam_path)]
                }
        } else {
            aligned_ch = Channel.empty()
        }
    }

    // ========================================================================
    // QC: Metrics collection and reporting
    // ========================================================================

    if (!params.skip_qc) {
        // Interop plot (run-level)
        INTEROP_PLOT(run_dir)

        // Per-sample QC metrics
        if (params.sampleinfo) {
            // Use sampleinfo for QC (required for panel-specific thresholds)
            qc_samples_ch = Channel
                .fromPath(params.sampleinfo)
                .splitCsv(sep: '\t', header: true)
                .map { row ->
                    def sample_id = row['Samples']
                    def panel = row['Panel'] ?: row['Tumor']
                    def bam_path = row['BAM Path']
                    [sample_id, panel, file(bam_path)]
                }

            PICARD_HS(qc_samples_ch)
            QCSUM_METRICS(PICARD_HS.out.hsmetrics)

            // Collect all qcsum files for MultiQC
            qcsum_ch = QCSUM_METRICS.out.qcsum.collect()
        } else {
            qcsum_ch = Channel.empty()
        }

        // MultiQC aggregation
        fastqs_dir = file("${params.fastqs_dir}/${run_name}")
        bams_dir = file(final_outdir)

        MULTIQC(
            fastqs_dir,
            bams_dir,
            qcsum_ch.ifEmpty([]),
            INTEROP_PLOT.out.plot.ifEmpty([])
        )
    }
}

// ----------------------------------------------------------------------------
// Completion handler
// ----------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    ════════════════════════════════════════════════════════════════
    Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration          : ${workflow.duration}
    Output directory  : ${final_outdir}
    ════════════════════════════════════════════════════════════════
    """.stripIndent()
}
