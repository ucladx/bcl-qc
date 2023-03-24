#!/bin/bash
set -eo pipefail

run_path=$1
samplesheet=$2
fastq_output=$3

if [ -d $fastq_output ]; then
  echo "Using existing folder for this runs's FASTQs: $fastq_output"
else
  echo "Creating new folder for this run's FASTQs: $fastq_output"
  mkdir $fastq_output
fi

echo "running demux on $samplesheet
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory $run_path --sample-sheet $samplesheet \
--output-directory $fastq_output
