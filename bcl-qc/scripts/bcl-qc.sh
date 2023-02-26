#!/bin/bash
set -eo pipefail

index=$1
run_name=$2
bed=$3

if [ $# -ne 3 ]; then
  echo "Usage: bash bcl-qc.sh <index> <run_name> <exon_bed>"
  exit 1
fi

echo "running bcl-qc.sh with index: $index"
### demux

echo "creating reads folder on dragen local disk"
mkdir /staging/hot/reads/$run_name

echo "running demux"
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory /mnt/pns/runs/$run_name --sample-sheet /mnt/pns/runs/$run_name/SampleSheet_$index.csv \
--output-directory /staging/hot/reads/$run_name/$index

### align

echo "creating bam dirs"
cut -f2 -d, /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} mkdir -p /mnt/pns/bams/$run_name/{}

echo "running alignment"
cut -f2 -d, /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} \
dragen --enable-map-align true --enable-map-align-output true --output-format BAM \
--enable-duplicate-marking true --generate-sa-tags true --enable-sort true --soft-read-trimmers polyg,quality --trim-min-quality 2 \
--ref-dir /staging/human/reference/hg38_alt_masked_graph_v2 --intermediate-results-dir /staging/tmp  \
--qc-coverage-tag-1 target_bed --qc-coverage-region-1 $bed  --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true \
--enable-variant-caller true --vc-combine-phased-variants-distance 6 --vc-emit-ref-confidence GVCF --enable-hla true \
--fastq-list /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv --fastq-list-sample-id {} \
--output-directory /mnt/pns/bams/$run_name/{} --output-file-prefix {}
