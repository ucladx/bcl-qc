#!/bin/bash
set -eo pipefail

run_path=$1
idx=$2
output_dir=$3
samplesheet=$run_path/SampleSheet_$idx.csv
fastq_output=$output_dir/$idx

mkdir -p $output_dir
echo "running demux using samplesheet: $samplesheet and fastq_output: $fastq_output"
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory $run_path --sample-sheet $samplesheet \
--output-directory $fastq_output
