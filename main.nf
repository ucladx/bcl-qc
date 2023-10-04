#!/usr/bin/env nextflow
include { validateParameters; paramsHelp; fromSamplesheet; paramsSummaryLog } from 'plugin/nf-validation'

process demux {
  input:
  path run_dir
  path samplesheet
  path fastq_outdir
  output:
  path fastq_list
  // Will dragen still output to e.g. /staging/hot/reads/$RUN_NAME/I10/... when we only use one samplesheet?
  // Otherwise need to parse samplesheet filename to be able to output the correct fastq_list path
  """
  dragen --bcl-conversion-only true \
  --bcl-use-hw false \
  --bcl-only-matched-reads true \
  --bcl-input-directory ${run_dir} \
  --sample-sheet ${samplesheet} \
  --output-directory ${fastq_outdir}

  echo ${fastq_outdir}$(basename ${run_dir})/Reports/fastq_list.csv
  """
}

process align {
  input:
  path fastq_list
  path bam_outdir
  tuple val(sample_id), val(assay)

  // assign variables dragen_install, reference_dir, targets_bed, baits_bed using config/bclqc.yaml and assay name
  script:
  if( assay == 'pcp' ) {

  }
  else if ( assy == 'pto') {

  }
  else if( assay == 'ces' ) {

  }
  else
    error "Invalid assay provided: ${assay}"
  // check dragen version (maybe just install ${dragen_install} every time?)
  // need sudo for dragen install? e.g.
  // sudo sh ${dragen_install}
  """
  mkdir ${bam_outdir}${sample_id}
  dragen \
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
}

process multiqc {
  input:
  path fastqs_dir
  path bams_dir
  
  """
  multiqc --force \
  --config config/multiqc_config.yaml \
  --outdir ${bams_dir} \
  ${bams_dir} \
  ${fastqs_dir}
  """
}

workflow {
  def samplesheet_info = Channel.fromSampleSheet("input").map { sample_id, barcode1, barcode2, assay -> tuple(sample_id, assay) }
  demux(params.run_dir, params.input) | align(params.bam_outdir, samplesheet_info) | multiqc(params.fastq_outdir, params.bam_outdir) | view
}

if (params.validate_params) {
  validateParameters()
}

if (params.show_params_summary) {
  log.info paramsSummaryLog(workflow)
}

WorkflowMain.initialise(workflow, params, log)
