index=$1
run_name=$2
EXONS=/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed

echo "creating bam dirs"
cut -f2 -d, /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} mkdir -p /mnt/pns/bams/$run_name/{}

echo "running alignment"
cut -f2 -d, /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} \
dragen --enable-map-align true --enable-map-align-output true --output-format BAM \
--enable-duplicate-marking true --generate-sa-tags true --enable-sort true \
--ref-dir /staging/human/reference/hg38_alt_masked_graph_v2 --intermediate-results-dir /staging/tmp  \
--qc-coverage-tag-1 target_bed --qc-coverage-region-1 $EXONS  --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true \
--enable-variant-caller true --vc-combine-phased-variants-distance 6 --vc-emit-ref-confidence GVCF --enable-hla true \
--fastq-list /staging/hot/reads/$run_name/$index/Reports/fastq_list.csv --fastq-list-sample-id {} \
--output-directory /mnt/pns/bams/$run_name/{} --output-file-prefix {}