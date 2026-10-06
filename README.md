# bcl-qc-nf

Nextflow pipeline for DRAGEN demultiplexing, alignment, Picard HsMetrics and MultiQC.
Outputs retain the run/sample layout of the original Python pipeline.

## Setup

On the DRAGEN server, install Conda in an executable location (such as `/opt/conda`),
then create and activate the pinned environment:

```bash
conda env create -f config/nf_env.yaml
conda activate bclqc-nf
```

`config/nf_env.yaml` is the sole dependency specification. The launcher checks its
pins against the activated environment. DRAGEN is supplied by the server. Adjust
reference paths in `nextflow.config`, QC settings in `config/qcsum_config.yaml`, and
report settings in `config/multiqc_config.yaml` before use.

## Production

Run from a release checkout inside tmux/screen:

```bash
scripts/run_bclqc_nf.sh --run_dir /path/to/RUN --sampleinfo /path/to/sampleinfo.tsv
```

The default steps are `demux,align,qc`; use `--steps align,qc` or `--steps qc` to rerun
selected steps. QC requires a TSV with unique `Samples`, `Panel` (or `Tumor`) and
`BAM Path` values, and an existing FASTQ run directory. Sample IDs may contain
letters, digits, underscores, dots and hyphens, starting with a letter or digit.
Run name is the basename of `--run_dir`.

| Parameter | Default |
|---|---|
| `--fastqs_dir` | `/staging/hot/reads` |
| `--bams_dir` | `/mnt/pns/bams` |
| `--dragen_tmp` | `/staging/tmp` |
| `--dragen_lock` | `/staging/tmp/.bclqc_dragen.lock` |
| `--dragen_lock_wait` | `43200` seconds |

Alignment handles `I10` (PCP) and `U5N2_I10` (HEME); `I8N2_N10` alignment is skipped.
Unknown indexes and duplicate samples across alignment indexes are refused.
Launches use a per-run lock; DRAGEN tasks share a host lock. Coordinate manual
DRAGEN commands separately because they do not acquire that lock.

FASTQs are written to `<fastqs_dir>/<RUN>/<index>/`; alignments and sample QC go to
`<bams_dir>/<RUN>/<sample>/`. Run outputs include `qcsum_mqc.csv`,
`occ_pf_lane_mqc.jpg` and `multiqc_report.html`. Occupancy plotting is optional.
Logs, versions, commit, trace, report and exit status go to `pipeline_info/` within
the run's BAM directory. Work uses `/staging/bclqc-nf-work` (`BCLQC_NF_WORK` to override)
and is removed on success; `BCLQC_NF_KEEP_WORK=1` retains it. Failed work is retained.

## Validation

From the activated environment, run scratch-only QC against existing production inputs:

```bash
scripts/dragen_validate.sh --run-dir /path/to/RUN --sampleinfo /path/to/sampleinfo.tsv
```

Results go to `/staging/bclqc-nf-validate/<RUN>/<timestamp>/bams/<RUN>/`. Override
paths with `--scratch`, `--prod-bams-run` and `--prod-fastq-run`. Scratch must be new
and separate from inputs and the repository. Existing sample files are linked into
scratch; Picard and qcsum are recomputed. Production inputs are read only.

Exit 0 means QC completed. Inspect the scratch `qcsum_mqc.csv` and `multiqc_report.html`
and compare with the production results before deployment. This validates the QC
path; it does not run DRAGEN demultiplexing/alignment or certify output equivalence.

For an end-to-end command review against a completed production run:

```bash
scripts/dragen_validate.sh --dry-run --run-dir /path/to/RUN --sampleinfo /path/to/sampleinfo.tsv
```

This prints and saves `commands.txt` under the fresh scratch directory. It lists
all demux, alignment, Picard, qcsum, merge and MultiQC commands as a full rerun,
even when production outputs already exist. Commands use the production paths
(selected with `--prod-bams-run` and `--prod-fastq-run`) for comparison with production
logs. Picard intervals and Perl QC thresholds are expanded. Exome alignment skips
are noted. The same command builders generate production task scripts and this report.

Dry-run reads the existing `Reports/fastq_list.csv` files for sample discovery and
requires complete demultiplexing for every sample sheet. It creates only scratch
metadata and the report: no DRAGEN, Picard, QC, deletion, or lock commands are executed.
The report assumes fresh outputs; actual runs still apply their existing-output
skip rules. This checks command construction, not DRAGEN execution or output quality.

## Recovery

Task caching is disabled; do not use `-resume`.

- **Demux:** remove each failed index directory and its adjacent
  `.<index>.bclqc_demux_incomplete` marker, then rerun demux. Complete indexes are
  skipped when others are missing. Alignment requires every expected FASTQ list.
- **Alignment:** remove directories containing `.bclqc_align_incomplete`, then rerun
  `--steps align,qc`. Existing alignment directories without a CRAM are refused.
  QC refuses incomplete alignments.
- **QC:** fix the cause and rerun `--steps qc`. Existing HsMetrics files are reused;
  delete them to recompute Picard. qcsum and MultiQC always rerun. Invalid or missing
  metrics fail explicitly. Picard writes a temporary file, then renames it on success.

Every QC run deletes `<bams_dir>/<RUN>/*/*.wgs_*.csv` before MultiQC, including on
historical runs. Deletions are logged in `pipeline_info/wgs_csv_deleted.txt`, with
no backup. Validation removes only scratch links. `--forbid_prefixes` protects
comma-separated paths from output writes, including through symlinks.
