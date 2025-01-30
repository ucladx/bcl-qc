#!/bin/bash

cram=$1
sample=$2

#conda activate bio

mkdir -p /mnt/pns/qc/heme/$sample

picard CollectHsMetrics \
I=$cram \
O=/mnt/pns/qc/heme/$sample.hsm.txt \
R=/mnt/pns/tracks/ref/hg38.fa \
BAIT_INTERVALS=/mnt/pns/tracks/goal_ucla_heme_221_merged_baits.hg38.ilist \
TARGET_INTERVALS=/mnt/pns/tracks/hemev2_roi.interval_list

perl qcsum_metrics.pl \
--prefix $sample \
--qcfolder "/mnt/pns/qc/heme/$sample" \
--pipeline_version v1.0d \
--platform NovaSeq6000 \
--pass_min_align_pct 99.0 \
--fail_min_align_pct 90.0 \
--covered 500 \
--pass_min_roi_pct 86 \
--fail_min_roi_pct 77 \
--pass_min_avgcov 1400 \
--fail_min_avgcov 760 \
--pass_min_reads 21000000 \
--fail_min_reads 12000000 \
--capture goal \
--capture_version v1.0
