#!/usr/bin/env bash
# Run scratch QC, or print a full demux/align/QC command plan without running analysis.
# Usage: scripts/dragen_validate.sh --run-dir DIR --sampleinfo TSV [options]
#   --scratch DIR         fresh output directory (default /staging/bclqc-nf-validate/RUN/TIMESTAMP)
#   --prod-bams-run DIR   existing alignment output (default /mnt/pns/bams/RUN)
#   --prod-fastq-run DIR  existing FASTQ output (default /staging/hot/reads/RUN)
#   --dry-run             print a full rerun using production paths and existing FASTQ lists
# Activate the environment from config/nf_env.yaml before running.
# Runs the current checkout; exit 0 means QC completed or the dry-run report was generated.
set -euo pipefail
HERE=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
RUN_DIR=""; SI=""; SCRATCH=""; PROD_BAMS=""; PROD_FQ=""; DRY_RUN=0
die() { echo "ERROR: $*" >&2; exit 2; }
inside() { [ "$2" = / ] || [ "$1" = "$2" ] || [[ $1 == "$2"/* ]]; }
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --run-dir|--sampleinfo|--scratch|--prod-bams-run|--prod-fastq-run)
      [ $# -ge 2 ] && [ -n "$2" ] || die "$1 needs a value"
      case "$1" in
        --run-dir) RUN_DIR=$2 ;;
        --sampleinfo) SI=$2 ;;
        --scratch) SCRATCH=$2 ;;
        --prod-bams-run) PROD_BAMS=$2 ;;
        --prod-fastq-run) PROD_FQ=$2 ;;
      esac
      shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done
[ -n "$RUN_DIR" ] && [ -d "$RUN_DIR" ] || die "--run-dir must name an existing directory"
[ -n "$SI" ] && [ -f "$SI" ] || die "--sampleinfo must name an existing TSV"
RUN_DIR=$(realpath -m -s "$RUN_DIR"); SI=$(readlink -f "$SI"); RUN=$(basename "$RUN_DIR")
SCRATCH=$(readlink -m "${SCRATCH:-/staging/bclqc-nf-validate/$RUN/$(date +%Y%m%d-%H%M%S)}")
PROD_BAMS=$(readlink -m "${PROD_BAMS:-/mnt/pns/bams/$RUN}")
PROD_FQ=$(realpath -m -s "${PROD_FQ:-/staging/hot/reads/$RUN}")
[ "$DRY_RUN" = 1 ] || [ -d "$PROD_BAMS" ] || die "alignment output not found: $PROD_BAMS"
[ "$DRY_RUN" = 0 ] || [ "$(basename "$PROD_BAMS")" = "$RUN" ] || die "--prod-bams-run must be named $RUN in dry-run mode"
[ -d "$PROD_FQ" ] && [ "$(basename "$PROD_FQ")" = "$RUN" ] || die "FASTQ run directory must exist and be named $RUN"
[ ! -e "$SCRATCH" ] && [ ! -L "$SCRATCH" ] || die "scratch directory already exists: $SCRATCH"
FORBID="/mnt/pns,/staging/hot,$PROD_BAMS,$PROD_FQ,$RUN_DIR"
for p in /mnt/pns /staging/hot "$PROD_BAMS" "$PROD_FQ" "$RUN_DIR" "$HERE" "$SI"; do
  p=$(readlink -m "$p")
  if inside "$SCRATCH" "$p" || inside "$p" "$SCRATCH"; then
    die "scratch directory overlaps an input or repository directory: $p"
  fi
done
[ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python3" ] || die "activate the environment from config/nf_env.yaml first"
mkdir -p "$(dirname "$SCRATCH")"
mkdir "$SCRATCH"
# Link only per-sample inputs. QC regenerates HsMetrics and qcsum; run-level outputs stay in scratch.
# Resolve relative BAM Paths before changing the launch directory.
"$CONDA_PREFIX/bin/python3" - "$SI" "$PROD_BAMS" "$SCRATCH/bams/$RUN" "$SCRATCH/sampleinfo.tsv" "$DRY_RUN" <<'PY'
import csv
import re
import sys
from pathlib import Path

src, prod, out, dest = map(Path, sys.argv[1:5])
dry_run = sys.argv[5] == "1"
with src.open(newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    fields = reader.fieldnames
    if not fields or not {'Samples', 'BAM Path'} <= set(fields) or not ({'Panel', 'Tumor'} & set(fields)):
        sys.exit('ERROR: sampleinfo needs Samples, Panel (or Tumor), and BAM Path columns')
    rows = list(reader)
if not rows:
    sys.exit('ERROR: sampleinfo has no samples')
for row in rows:
    sample = row['Samples']
    if not sample or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample) or sample == 'pipeline_info' or sample.startswith('multiqc'):
        sys.exit(f'ERROR: invalid sample directory name: {sample!r}')
    if not row['BAM Path']:
        sys.exit(f'ERROR: missing BAM Path for {sample}')
    bam = Path(row['BAM Path']).absolute()
    if not dry_run and (not bam.is_file() or (bam.parent / '.bclqc_align_incomplete').exists()):
        sys.exit(f'ERROR: missing or incomplete alignment: {bam}')
    row['BAM Path'] = str(bam)
samples = {row['Samples'] for row in rows}
if len(samples) != len(rows):
    sys.exit('ERROR: duplicate samples in sampleinfo')
if not dry_run:
    for directory in prod.iterdir():
        sample = directory.name
        if not directory.is_dir() or sample == 'pipeline_info' or sample.startswith('multiqc'):
            continue
        if sample not in samples and not (directory / f'{sample}.cram').is_file():
            continue
        if (directory / '.bclqc_align_incomplete').exists():
            sys.exit(f'ERROR: incomplete alignment: {directory}')
        target = out / sample
        target.mkdir(parents=True)
        for path in directory.iterdir():
            if path.is_file() and not path.name.endswith(('.qcsum.txt', '.hsm.txt', '.hsm.txt.partial')):
                (target / path.name).symlink_to(path)
with dest.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
PY
cd "$SCRATCH"
if [ "$DRY_RUN" = 1 ]; then
  export PATH="$CONDA_PREFIX/bin:$PATH" JAVA_HOME="$CONDA_PREFIX/lib/jvm" PYTHONNOUSERSITE=1 NXF_ANSI_LOG=false
  unset JAVA_CMD NXF_JAVA_HOME PYTHONPATH PYTHONHOME
  exec nextflow -log "$SCRATCH/nextflow.log" run "$HERE" -work-dir "$SCRATCH/work" \
    --dry_run true --steps demux,align,qc --run_dir "$RUN_DIR" --sampleinfo "$SCRATCH/sampleinfo.tsv" \
    --fastqs_dir "$(dirname "$PROD_FQ")" --bams_dir "$(dirname "$PROD_BAMS")"
fi
export BCLQC_NF_WORK="$SCRATCH/work"
echo "QC results: $SCRATCH/bams/$RUN"
exec "$HERE/scripts/run_bclqc_nf.sh" --run_dir "$RUN_DIR" --steps qc \
  --sampleinfo "$SCRATCH/sampleinfo.tsv" --fastqs_dir "$(dirname "$PROD_FQ")" \
  --bams_dir "$SCRATCH/bams" --dragen_tmp "$SCRATCH/dragen_tmp" --forbid_prefixes "$FORBID"
