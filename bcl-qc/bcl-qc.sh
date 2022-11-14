run_dir=$1

# demux
echo "copying $run_dir to staging"
# rsync -rhlW --info=progress2 --exclude='Thumbnail_Images' /mnt/pns/runs/$run_dir /staging/hot/$run_dir

mkdir /staging/hot/reads/$run_dir
# mkdir /staging/hot/reads/$run_dir/I10

dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true --bcl-input-directory /staging/hot/$run_dir --sample-sheet /staging/hot/$run_dir/SampleSheet_I10.csv --output-directory /staging/hot/reads/$run_dir/I10

# align
cut -f2 -d, /staging/hot/reads/$run_dir/I10/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} mkdir -p /mnt/pns/bams/$run_dir/{}
cut -f2 -d, /staging/hot/reads/$run_dir/I10/Reports/fastq_list.csv | grep -v ^RGSM | sort -u | xargs -L1 -I{} dragen --intermediate-results-dir /staging/tmp --enable-map-align true --enable-map-align-output true --output-format BAM --enable-duplicate-marking true --generate-sa-tags true --enable-sort true --ref-dir /staging/human/reference/hg38_alt_masked_graph_v2 --qc-coverage-tag-1 target_bed --qc-coverage-region-1 $EXONS --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true --enable-variant-caller true --vc-combine-phased-variants-distance 6 --vc-emit-ref-confidence GVCF --enable-hla true --fastq-list /staging/hot/reads/$run_dir/I10/Reports/fastq_list.csv --fastq-list-sample-id {} --output-directory /mnt/pns/bams/$run_dir/{} --output-file-prefix {}

# multiqc
rm -f /mnt/pns/bams/$run_dir/*/*.wgs_*.csv
multiqc --outdir /mnt/pns/bams/$run_dir /mnt/pns/bams/$run_dir /staging/hot/reads/$run_dir/I10