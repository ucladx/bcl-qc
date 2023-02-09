index=$1
run_name=$2
echo "running multiqc.sh"

rm -f /mnt/pns/bams/$run_name/*/*.wgs_*.csv
multiqc --force --config ../config/multiqc_config.yaml --outdir ~ /mnt/pns/bams/$run_name /staging/hot/reads/$run_name/$index/
