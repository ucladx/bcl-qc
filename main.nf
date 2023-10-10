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
  each sample_id

  when:
  sample_id != null

  script:
  // TODO validate assay names using samplesheet schema
  sample_id = sample_id.trim()
  assay_params = params[params.assay]
  dragen_install = assay_params.dragen_install
  reference_dir = assay_params.reference_dir
  targets_bed = assay_params.targets_bed
  baits_bed = assay_params.baits_bed
  addtl_cmds = assay_params.addtl_cmds

  outdir = "${params.bam_outdir}" + "${sample_id}"

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
--output-directory ${outdir} \
--output-file-prefix ${sample_id}"""

  if (addtl_cmds) {
    alignment_cmd += " ${addtl_cmds}"
  }

  shell:
  '''
  mkdir -p !{outdir}
  !{alignment_cmd}
  '''
}

process multiqc {
  input:
  path fastqs_dir
  path bams_dir
  
  shell:
  '''multiqc --force \
  --config config/multiqc_config.yaml \
  --outdir !{params.bams_dir} \
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
  start_line=$(grep -n 'BCLConvert_Data' "!{samplesheet}" | cut -d: -f1)
  if [[ -n "$start_line" ]]; then
      sed -n "$((start_line + 1)),$ p" "!{samplesheet}" > 'trimmed_samplesheet.csv'
  fi
  '''
}

process parse_fastq_list {
  input:
  path fastq_list

  output:
  path 'samples.txt'

  shell:
  '''
  cat !{fastq_list} | cut -d, -f2 | tail +2 > samples.txt
  '''
}

workflow {
  fastq_list = params.fastq_list
  log.info paramsSummaryLog(workflow)

  if ('demux' in params.steps) {
    fastq_list = demux(params.run_dir, params.input, params.fastq_outdir)
  }

  if ('align' in params.steps) {
    if (params.assay) {
      if (fastq_list) {
        align(fastq_list, params.bam_outdir, parse_fastq_list(fastq_list).splitText())
      }
    }
    else {
      def samplesheet_info = trim_samplesheet(params.input).map{sample,index,index2,assay -> tuple(sample,params.assay)}
      align(fastq_list, params.bam_outdir, samplesheet_info)
    }
  }

  if ('multiqc' in params.steps) {
    multiqc(params.fastq_outdir, params.bam_outdir)
  }
}
