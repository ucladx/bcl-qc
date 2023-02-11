index=$1
run_name=$2

echo "creating reads folder on dragen local disk"
mkdir /staging/hot/reads/$run_name

echo "running demux"
dragen --bcl-conversion-only true --bcl-use-hw false --bcl-only-matched-reads true \
--bcl-input-directory /mnt/pns/runs/$run_name --sample-sheet /mnt/pns/runs/$run_name/SampleSheet_$index.csv \
--output-directory /staging/hot/reads/$run_name/$index