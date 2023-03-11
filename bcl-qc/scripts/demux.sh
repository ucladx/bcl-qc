#!/bin/bash
set -eo pipefail

index=$1
run_name=$2

if [ $# -ne 2 ]; then
  echo "Usage: bash demux.sh <index> <run_name>"
  exit 1
fi

fastq_dir="/staging/hot/reads/$run_name"
if [ -d $fastq_dir ]; then
  echo "Using existing folder for this runs's FASTQs: $fastq_dir"
else
  echo "Creating new folder for this run's FASTQs: $fastq_dir"
  mkdir $fastq_dir
fi

echo "running demux"
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory /mnt/pns/runs/$run_name --sample-sheet /mnt/pns/runs/$run_name/SampleSheet_$index.csv \
--output-directory /staging/hot/reads/$run_name/$index
