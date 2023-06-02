#!/bin/bash
set -eo pipefail

fastq_list=$1
output=$2
bed=$3
sample_id=$4

mkdir -p $output
echo "running align on $sample_id"
dragen --enable-map-align true --enable-map-align-output true --output-format BAM \
--enable-duplicate-marking true --generate-sa-tags true --enable-sort true --soft-read-trimmers polyg,quality --trim-min-quality 2 \
--ref-dir /staging/human/reference/hg38_alt_masked_graph_v2 --intermediate-results-dir /staging/tmp  \
--qc-coverage-tag-1 target_bed --qc-coverage-region-1 $bed  --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true \
--enable-variant-caller true --vc-combine-phased-variants-distance 6 --vc-emit-ref-confidence GVCF --enable-hla true \
--fastq-list $fastq_list --fastq-list-sample-id $sample_id \
--output-directory $output --output-file-prefix $sample_id
