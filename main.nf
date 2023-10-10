#!/usr/bin/env nextflow
include { validateParameters; paramsHelp; fromSamplesheet; paramsSummaryLog } from 'plugin/nf-validation'

process demux {
  input:
  path run_name
  path samplesheet
  path fastq_outdir
  output:
  path "fastq_list.csv"

  shell:
  '''
  dragen --bcl-conversion-only true \
  --bcl-use-hw false \
  --bcl-only-matched-reads true \
  --bcl-input-directory !{params.run_dir} \
  --sample-sheet !{samplesheet} \
  --output-directory !{params.fastq_outdir}

  cat !{params.fastq_outdir}!{run_name}/Reports/fastq_list.csv > fastq_list.csv
  '''
}

process align {
  input:
  path fastq_list
  path bam_outdir
  tuple val(sample_id), val(assay)

  script:
  // TODO validate assay names using samplesheet schema
  assay_params = params[assay]
  dragen_install = assay_params.dragen_install
  reference_dir = assay_params.reference_dir
  targets_bed = assay_params.targets_bed
  baits_bed = assay_params.baits_bed

  alignment_cmd = """dragen \
  --enable-map-align true \
  --enable-map-align-output true \
  --output-format BAM \
  --enable-duplicate-marking true \
  --generate-sa-tags true \
  --enable-sort true \
  --soft-read-trimmers polyg,quality \
  --trim-min-quality 2 \
  --ref-dir ${reference_dir} \
  --intermediate-results-dir /staging/tmp \
  --qc-coverage-tag-1 target_bed \
  --qc-coverage-region-1 ${targets_bed} \
  --qc-coverage-reports-1 cov_report \
  --qc-coverage-ignore-overlaps true \
  --enable-variant-caller true \
  --vc-combine-phased-variants-distance 6 \
  --vc-emit-ref-confidence GVCF \
  --enable-hla true \
  --fastq-list ${fastq_list} \
  --fastq-list-sample-id ${sample_id} \
  --output-directory ${bam_outdir} \
  --output-file-prefix ${sample_id}
  """

  if (assay_params.dragen_addtl_align_opts) {
    alignment_cmd += assay_params.dragen_addtl_align_opts
  }

  shell:
  '''
  sh !{dragen_install}
  mkdir !{bam_outdir}!{sample_id}
  !{alignment_cmd}
  '''
}

process multiqc {
  input:
  path fastqs_dir
  path bams_dir
  
  shell:
  '''
  multiqc --force \
  --config config/multiqc_config.yaml \
  --outdir !{bams_dir} \
  !{bams_dir} \
  !{fastqs_dir}
  '''
}

process trim_samplesheet{
  input:
  path samplesheet

  output:
  path 'trimmed_samplesheet.csv'

  shell:
  '''
  awk '[BCLConvert_Data] {flag=1; next} flag' !{samplesheet} > trimmed_samplesheet.csv
  '''
}

workflow {
  def fastq_list = params.fastq_list // null if not provided
  def dragen_samplesheet = params.input
  if ('demux' in params.steps) {
    demux(params.run_dir, dragen_samplesheet, params.fastq_outdir)
  }
  if ("align" in params.steps) {
    def samplesheet_info = Channel.fromSamplesheet(trim_samplesheet(dragen_samplesheet)).map{sample,index,index2,assay -> tuple(sample,assay)}
    align(fastq_list, params.bam_outdir, samplesheet_info)
  }
  if ("multiqc" in params.steps) {
    multiqc(params.fastq_outdir, params.bam_outdir)
  }
  log.info paramsSummaryLog(workflow)
}

// if (params.validate_params) {
//   validateParameters()
// }
