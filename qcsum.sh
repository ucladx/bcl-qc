#!/bin/bash

cram=$1
<<<<<<< HEAD
sample=$2
out_dir=$3

mkdir -p $out_dir

picard CollectHsMetrics \
I=$cram \
O=$out_dir/$sample.hsm.txt \
R=/mnt/pns/tracks/ref/hg38.fa \
BAIT_INTERVALS=/mnt/pns/tracks/goal_ucla_heme_221_merged_baits.hg38.ilist \
TARGET_INTERVALS=/mnt/pns/tracks/hemev2_roi.interval_list

perl qcsum_metrics.pl \
--prefix $sample \
--qcfolder "$out_dir" \
--pipeline_version HemeDx-v2.0 \
=======
output=$2

#conda activate bio

picard CollectHsMetrics \
I=$cram \
O=/mnt/pns/qc/heme/$output.hsm.txt \
R=/mnt/pns/ref/hg38.fa \
BAIT_INTERVALS=/mnt/pns/tracks/goal_ucla_heme_221_merged_baits.hg38.ilist \
TARGET_INTERVALS=/mnt/pns/tracks/hemev2_roi.interval_list

sample=$output
perl qcsum_metrics.pl \
--prefix $sample \
--qcfolder "/mnt/pns/qc/heme/" \
--pipeline_version v1.0d \
>>>>>>> b92b6e8 (Add qcsum scripts)
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
<<<<<<< HEAD
=======

>>>>>>> b92b6e8 (Add qcsum scripts)
