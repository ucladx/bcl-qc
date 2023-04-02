#!/bin/bash
set -eo pipefail

run_path=$1
samplesheet=$2
fastq_output=$3

echo "running demux on $samplesheet
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory $run_path --sample-sheet $samplesheet \
--output-directory $fastq_output
