#!/usr/bin/env bash
# Usage: conda activate bclqc-nf; scripts/run_bclqc_nf.sh --run_dir DIR --sampleinfo TSV [pipeline options]
# BCLQC_NF_WORK sets the work base (default /staging/bclqc-nf-work); BCLQC_NF_KEEP_WORK=1 retains successful work.
# Run inside tmux/screen. Logs and exit status are saved under <bams_dir>/<RUN>/pipeline_info/.
set -euo pipefail
HERE=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
RUN_DIR=""; BAMS_DIR=/mnt/pns/bams; FORBID=""; STEPS=demux,align,qc
ARGS=("$@")
die() { echo "ERROR: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help) sed -n '2,4p' "$0"; exit 0 ;;
    --dry_run|--dry_run=*) die "use scripts/dragen_validate.sh --dry-run for command planning" ;;
    --run_dir|--bams_dir|--forbid_prefixes|--steps)
      [ $# -ge 2 ] && [[ $2 != --* ]] || die "$1 needs a value"
      key=$1; value=$2; shift ;;
    --run_dir=*|--bams_dir=*|--forbid_prefixes=*|--steps=*) key=${1%%=*}; value=${1#*=} ;;
    *) shift; continue ;;
  esac
  case "$key" in
    --run_dir) RUN_DIR=$value ;; --bams_dir) BAMS_DIR=$value ;;
    --forbid_prefixes) FORBID=$value ;; --steps) STEPS=$value ;;
  esac
  shift
done
[ -d "$RUN_DIR" ] || die "--run_dir must name an existing directory"
[ -n "$BAMS_DIR" ] || die "--bams_dir must not be empty"
[ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python3" ] || die "activate the environment from config/nf_env.yaml first"
export PATH="$CONDA_PREFIX/bin:$PATH" JAVA_HOME="$CONDA_PREFIX/lib/jvm" NXF_ANSI_LOG=false
unset JAVA_CMD NXF_JAVA_HOME PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
for tool in nextflow java picard multiqc perl python3; do
  [ "$(command -v "$tool" || true)" = "$CONDA_PREFIX/bin/$tool" ] || die "$tool is missing from $CONDA_PREFIX"
done
command -v flock >/dev/null || die "flock is required"
if [[ $STEPS == *demux* || $STEPS == *align* ]]; then
  command -v dragen >/dev/null || die "dragen is required for demux/alignment"
fi
# Check the activated environment against the same pins used to create it.
VERSIONS=$(python3 - "$HERE/config/nf_env.yaml" <<'PY'
import json, os, sys
from importlib.metadata import version
from pathlib import Path
import yaml

packages = {}
for path in (Path(os.environ['CONDA_PREFIX']) / 'conda-meta').glob('*.json'):
    item = json.loads(path.read_text())
    packages[item['name']] = item['version']
for dependency in yaml.safe_load(Path(sys.argv[1]).read_text())['dependencies']:
    pins = dependency['pip'] if isinstance(dependency, dict) else [dependency]
    for pin in pins:
        name, sep, wanted = pin.replace('==', '=').partition('=')
        if not sep:
            continue
        actual = version(name) if isinstance(dependency, dict) else packages.get(name)
        if actual != wanted:
            sys.exit(f'ERROR: {name}: found {actual}, expected {wanted}; recreate the environment from config/nf_env.yaml')
        print(f'{name}: {actual}')
PY
)
RUN=$(basename "$(realpath -m -s "$RUN_DIR")")
OUT=$(readlink -m "$BAMS_DIR")/$RUN
PI=$OUT/pipeline_info
for path in "$PI" "$PI/.lock"; do
  [ ! -L "$path" ] || die "refusing symlink: $path"
done
python3 "$HERE/bin/check_write_paths.py" "$FORBID" "$OUT" "$PI" "$PI/.lock"
mkdir -p "$PI"
exec 9>>"$PI/.lock"
flock -n 9 || { echo "ERROR: another launch for $RUN holds $PI/.lock" >&2; exit 3; }
INFO=$(mktemp -d "$PI/$(date +%Y%m%d_%H%M%S).XXXXXX")
NFPID=""; WORK=""
finish() {
  local rc=$?
  trap - EXIT INT TERM HUP
  if [ -n "$NFPID" ]; then kill -TERM "$NFPID" 2>/dev/null || true; wait "$NFPID" 2>/dev/null || true; fi
  echo "exit: $rc   finished: $(date -Is)" >> "$INFO/versions.txt"
  if [ -n "$WORK" ] && [ "$rc" = 0 ] && [ "${BCLQC_NF_KEEP_WORK:-0}" != 1 ] && [ ! -L "$WORK" ] && [ "$(dirname "$WORK")" = "$WORK_BASE" ]; then
    rm -rf -- "$WORK"
  fi
  echo "bcl-qc-nf: exit $rc; logs: $INFO"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP
WORK_BASE=${BCLQC_NF_WORK:-/staging/bclqc-nf-work}
python3 "$HERE/bin/check_write_paths.py" "$FORBID" "$WORK_BASE"
mkdir -p "$WORK_BASE"
WORK_BASE=$(readlink -f "$WORK_BASE")
WORK=$(mktemp -d "$WORK_BASE/$RUN.XXXXXX")
{
  echo "date: $(date -Is)   host: $(hostname)   user: $(id -un)"
  echo "commit: $(git -C "$HERE" rev-parse HEAD)"
  git -C "$HERE" status --short
  printf 'command: %q ' "$0"; printf '%q ' "${ARGS[@]}"; echo
  echo "env: $CONDA_PREFIX   work: $WORK"
  echo "$VERSIONS"
  if command -v dragen >/dev/null; then dragen --version || true; fi
} > "$INFO/versions.txt" 2>&1

nextflow -log "$INFO/nextflow.log" run "$HERE" -work-dir "$WORK" \
  -with-trace "$INFO/trace.txt" -with-report "$INFO/report.html" \
  "${ARGS[@]}" > >(tee "$INFO/console.log") 2>&1 &
NFPID=$!
RC=0; wait "$NFPID" || RC=$?
NFPID=""
exit "$RC"
