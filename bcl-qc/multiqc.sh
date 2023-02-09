run_name=$1
echo "running multiqc.sh"

rm -f /mnt/pns/bams/$run_name/*/*.wgs_*.csv
multiqc --force --config ../config/multiqc_config.yaml --outdir ~/home/iatol/ /mnt/pns/bams/$run_name /staging/hot/reads/$run_name/I10/