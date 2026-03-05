#!/usr/bin/env python3
"""
This script provides a convenient CLI to run the picard/qcsum workflow on
independent samples (CRAM files) with panel-specific cutoffs.

Usage:
  # Run on a PCP tumor sample
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_tumor

  # Run on a Heme sample with custom sample ID
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_comp --sample-id MySample

  # Run on a Heme subpanel (MPN)
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_mpn

  # Output to a specific directory
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_normal --output-dir /path/to/output

  # List available panels
  python3 picard_qcsum.py --list-panels

  # Merge multiple qcsum files
  python3 picard_qcsum.py --merge sample1.qcsum.txt sample2.qcsum.txt -o /path/to/output

  # Batch process multiple CRAMs and merge results
  python3 picard_qcsum.py --batch sample1.cram sample2.cram --panel pcp_tumor -o /path/to/output
"""

import sys
import argparse
import subprocess
import logging
import yaml
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
QCSUM_CONFIG_YAML = "/home/iatol/bcl-qc/config/qcsum_config.yaml"
PICARD_REF = "/home/iatol/hg38.fa"

QC_SUM_HEADER = (
    "Sample,Sequencing_Platform,Pipeline_version,Alignment_QC,Coverage_QC,"
    "Total_Reads,%Reads_Aligned,Capture,Avg_Capture_Coverage,%On/Near_Bait_Bases,"
    "%On_Bait_Bases,%On_Target_Bases,FOLD_80_BASE_PENALTY,AT_DROPOUT,GC_DROPOUT,"
    "Avg_ROI_Coverage,MEDIAN_ROI_COVERAGE,MAX_ROI_COVERAGE,"
    "%ROI_1x,%ROI_20x,%ROI_100x,%ROI_250x,%ROI_500x"
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


class QCSumConfig:
    """Manages QC summary configuration for a given panel."""

    VALID_PANELS = [
        "pcp_tumor", "pcp_normal",
        "heme_comp", "heme_mpn", "heme_jak2", "heme_pblp"
    ]

    def __init__(self, panel: str, config_path: Path = QCSUM_CONFIG_YAML):
        self.panel = panel
        self.config_path = config_path
        self.config = self._load_config()

    def _load_config(self) -> dict:
        """Load and validate panel configuration from YAML."""
        with open(self.config_path) as f:
            yaml_obj = yaml.safe_load(f)
            if not yaml_obj:
                raise ValueError(f"Failed to load {self.config_path}")

            if self.panel not in yaml_obj:
                available = ", ".join(yaml_obj.keys())
                raise ValueError(
                    f"Panel '{self.panel}' not found. Available panels: {available}"
                )

            # Heme subpanels inherit from heme_comp with overridden target intervals
            if self.panel.startswith("heme_") and self.panel != "heme_comp":
                yaml_dict = yaml_obj.get("heme_comp", {}).copy()
                subpanel_dict = yaml_obj.get(self.panel, {})
                if "target_intervals" in subpanel_dict:
                    yaml_dict["target_intervals"] = subpanel_dict["target_intervals"]
            else:
                yaml_dict = yaml_obj.get(self.panel, {})

            return {k: str(v) for k, v in yaml_dict.items()}

    def get(self, key: str, default: str = "") -> str:
        """Get a config value by key."""
        return self.config.get(key, default)


def run_command(cmd: list, dry_run: bool = False) -> subprocess.CompletedProcess:
    """Execute a command with logging and error handling."""
    cmd_str = ' '.join(cmd)
    logging.info(f"Running: {cmd_str}")

    if dry_run:
        logging.info("[DRY RUN] Command not executed")
        return subprocess.CompletedProcess(cmd, 0, b"", b"")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.stderr:
            logging.warning(f"stderr: {result.stderr.decode()}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed (exit {e.returncode}): {cmd_str}")
        if e.stderr:
            logging.error(f"stderr: {e.stderr.decode()}")
        raise


def run_picard_hs_metrics(
    cram_file: str,
    output_file: str,
    bait_intervals: str,
    target_intervals: str,
    reference: str = PICARD_REF,
    dry_run: bool = False
) -> None:
    """Run Picard CollectHsMetrics on a CRAM file."""
    picard_cmd = [
        "picard", "CollectHsMetrics",
        f"I={cram_file}",
        f"O={output_file}",
        f"R={reference}",
        f"BAIT_INTERVALS={bait_intervals}",
        f"TARGET_INTERVALS={target_intervals}"
    ]
    run_command(picard_cmd, dry_run)


def parse_hs_metrics(hsm_file: str) -> dict:
    """
    Parse Picard HsMetrics output file.

    The metrics are on line 7 of the .hsm.txt file (0-indexed),
    with headers on line 6.
    """
    with open(hsm_file) as f:
        lines = f.readlines()

    # Find the METRICS section - headers followed by values
    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET"):
            headers = line.strip().split('\t')
            values = lines[i + 1].strip().split('\t')
            return dict(zip(headers, values))

    raise ValueError(f"Could not parse metrics from {hsm_file}")


def generate_qcsum_line(
    sample_id: str,
    metrics: dict,
    config: QCSumConfig
) -> str:
    """
    Generate a single line for qcsum output in CSV format.

    Calculates alignment and coverage QC PASS/FAIL status and builds
    the CSV line matching QC_SUM_HEADER.
    """
    total_reads = float(metrics.get("TOTAL_READS", 0))
    mean_target_coverage = float(metrics.get("MEAN_TARGET_COVERAGE", 0))
    pct_aligned = float(metrics.get("PCT_PF_UQ_READS_ALIGNED", 0)) * 100

    covered_threshold = int(config.get("covered", 100))
    roi_coverage_col = f"PCT_TARGET_BASES_{covered_threshold}X"
    pct_roi_at_threshold = float(metrics.get(roi_coverage_col, 0)) * 100

    # Determine alignment QC (PASS/FAIL)
    fail_min_align = float(config.get("fail_min_align_pct", 0))
    fail_min_reads = float(config.get("fail_min_reads", 0))

    if pct_aligned < fail_min_align or total_reads < fail_min_reads:
        alignment_qc = "FAIL"
    else:
        alignment_qc = "PASS"

    # Determine coverage QC (PASS/FAIL)
    fail_min_roi = float(config.get("fail_min_roi_pct", 0))
    fail_min_avgcov = float(config.get("fail_min_avgcov", 0))

    if pct_roi_at_threshold < fail_min_roi or mean_target_coverage < fail_min_avgcov:
        coverage_qc = "FAIL"
    else:
        coverage_qc = "PASS"

    def get_pct(key: str) -> str:
        """Convert Picard fraction to percentage string."""
        val = float(metrics.get(key, 0))
        return f"{val * 100:.2f}"

    fields = [
        sample_id,
        config.get("platform", "NovaSeq6000"),
        config.get("pipeline_version", "Unknown"),
        alignment_qc,
        coverage_qc,
        metrics.get("TOTAL_READS", "0"),
        f"{pct_aligned:.2f}",
        f"{config.get('capture', 'unknown')}-{config.get('capture_version', 'v1')}",
        f"{float(metrics.get('MEAN_BAIT_COVERAGE', 0)):.2f}",
        get_pct("PCT_SELECTED_BASES"),  # % on/near bait bases
        get_pct("PCT_USABLE_BASES_ON_BAIT"),  # % on bait bases
        get_pct("PCT_USABLE_BASES_ON_TARGET"),  # % on target bases
        f"{float(metrics.get('FOLD_80_BASE_PENALTY', 0)):.2f}",
        f"{float(metrics.get('AT_DROPOUT', 0)):.2f}",
        f"{float(metrics.get('GC_DROPOUT', 0)):.2f}",
        f"{mean_target_coverage:.2f}",
        metrics.get("MEDIAN_TARGET_COVERAGE", "0"),
        metrics.get("MAX_TARGET_COVERAGE", "0"),
        get_pct("PCT_TARGET_BASES_1X"),
        get_pct("PCT_TARGET_BASES_20X"),
        get_pct("PCT_TARGET_BASES_100X"),
        get_pct("PCT_TARGET_BASES_250X"),
        get_pct("PCT_TARGET_BASES_500X"),
    ]

    return ",".join(fields)


def run_qcsum(
    cram_file: str,
    panel: str,
    sample_id: str = None,
    output_dir: str = None,
    reference: str = PICARD_REF,
    dry_run: bool = False,
) -> dict:
    """
    Run the complete picard/qcsum workflow on a single sample.

    Args:
        cram_file: Path to input CRAM/BAM file.
        panel: Panel name for selecting cutoffs (e.g., 'pcp_tumor', 'heme_comp').
        sample_id: Sample identifier. If not provided, derived from filename.
        output_dir: Directory for output files. Defaults to same directory as CRAM.
        reference: Path to reference FASTA for Picard.
        dry_run: If True, print commands without executing.

    Returns:
        dict with QC results and output file paths.
    """
    cram_path = Path(cram_file).resolve()

    if not cram_path.exists() and not dry_run:
        raise FileNotFoundError(f"CRAM file not found: {cram_path}")

    # Derive sample ID from filename if not provided
    if sample_id is None:
        sample_id = cram_path.stem
        if sample_id.endswith(".cram") or sample_id.endswith(".bam"):
            sample_id = Path(sample_id).stem

    # Set output directory
    if output_dir is None:
        output_dir = cram_path.parent
    else:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Processing sample: {sample_id}")
    logging.info(f"Input CRAM: {cram_path}")
    logging.info(f"Panel: {panel}")
    logging.info(f"Output directory: {output_dir}")

    # Load panel configuration
    config = QCSumConfig(panel)

    bait_intervals = config.get("bait_intervals")
    target_intervals = config.get("target_intervals")

    if not bait_intervals or not target_intervals:
        raise ValueError(f"Panel '{panel}' missing bait_intervals or target_intervals")

    # Define output files
    hsm_file = output_dir / f"{sample_id}.hsm.txt"
    qcsum_file = output_dir / f"{sample_id}.qcsum.txt"

    # Run Picard CollectHsMetrics (skip if output already exists)
    if hsm_file.exists():
        logging.info(f"Picard output exists, skipping: {hsm_file}")
    else:
        run_picard_hs_metrics(
            cram_file=str(cram_path),
            output_file=str(hsm_file),
            bait_intervals=bait_intervals,
            target_intervals=target_intervals,
            reference=reference,
            dry_run=dry_run
        )

    if dry_run:
        logging.info("[DRY RUN] Skipping metrics parsing")
        return {"sample_id": sample_id, "dry_run": True}

    # Parse Picard metrics
    logging.info(f"Parsing metrics from: {hsm_file}")
    metrics = parse_hs_metrics(str(hsm_file))

    # Generate qcsum output
    qcsum_line = generate_qcsum_line(sample_id, metrics, config)

    # The qcsum line is CSV; alignment_qc and coverage_qc are fields 3 and 4
    qcsum_fields = qcsum_line.split(",")
    alignment_qc = qcsum_fields[3]
    coverage_qc = qcsum_fields[4]

    pct_aligned = float(metrics.get("PCT_PF_UQ_READS_ALIGNED", 0)) * 100
    total_reads = int(metrics.get("TOTAL_READS", 0))
    mean_target_coverage = float(metrics.get("MEAN_TARGET_COVERAGE", 0))
    covered_threshold = int(config.get("covered", 100))
    pct_roi = float(metrics.get(f"PCT_TARGET_BASES_{covered_threshold}X", 0)) * 100

    logging.info(f"Alignment QC: {alignment_qc}")
    logging.info(f"  Total reads: {total_reads:,}")
    logging.info(f"  % Aligned: {pct_aligned:.2f}%")
    logging.info(f"Coverage QC: {coverage_qc}")
    logging.info(f"  Mean target coverage: {mean_target_coverage:.2f}x")
    logging.info(f"  % ROI at {covered_threshold}x: {pct_roi:.2f}%")

    with open(qcsum_file, "w") as f:
        f.write(QC_SUM_HEADER + "\n")
        f.write(qcsum_line + "\n")

    logging.info(f"QC summary written to: {qcsum_file}")

    return {
        "sample_id": sample_id,
        "panel": panel,
        "hsm_file": str(hsm_file),
        "qcsum_file": str(qcsum_file),
        "alignment_qc": alignment_qc,
        "coverage_qc": coverage_qc,
    }


def merge_qcsum_files(
    qcsum_files: list,
    output_dir: str = None,
    output_name: str = "merged_qcsum.csv"
) -> str:
    """
    Merge multiple qcsum files into a single file.

    Args:
        qcsum_files: List of paths to qcsum files to merge.
        output_dir: Directory for output file. Defaults to current directory.
        output_name: Name of the merged output file.

    Returns:
        Path to the merged output file.
    """
    if not qcsum_files:
        raise ValueError("No qcsum files provided to merge")

    # Set output directory
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / output_name

    all_lines = []
    header = None

    for qcsum_file in qcsum_files:
        qcsum_path = Path(qcsum_file)
        if not qcsum_path.exists():
            logging.warning(f"Skipping missing file: {qcsum_path}")
            continue

        with open(qcsum_path) as f:
            lines = f.readlines()

        if not lines:
            logging.warning(f"Skipping empty file: {qcsum_path}")
            continue

        # First file sets the header
        if header is None:
            header = lines[0].strip()

        # Add data lines (skip header)
        for line in lines[1:]:
            line = line.strip()
            if line:
                all_lines.append(line)

    if header is None:
        raise ValueError("No valid qcsum files found to merge")

    # Write merged output
    with open(output_file, "w") as f:
        f.write(header + "\n")
        for line in all_lines:
            f.write(line + "\n")

    logging.info(f"Merged {len(all_lines)} samples from {len(qcsum_files)} files")
    logging.info(f"Output: {output_file}")

    return str(output_file)


def list_panels():
    """List available panels from the config file."""
    with open(QCSUM_CONFIG_YAML) as f:
        yaml_obj = yaml.safe_load(f)

    print("Available panels:")
    for panel in yaml_obj.keys():
        config = yaml_obj[panel]
        if "pipeline_version" in config:
            print(f"  {panel:15s} - {config.get('pipeline_version', '')}")
            print(f"    Coverage threshold: {config.get('covered', 'N/A')}x")
            print(f"    Min avg coverage: {config.get('pass_min_avgcov', 'N/A')}x")
        else:
            # Subpanel - inherits from parent
            print(f"  {panel:15s} - (subpanel, inherits from heme_comp)")
    return yaml_obj.keys()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Picard CollectHsMetrics and QC summary on a CRAM file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on a PCP tumor sample
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_tumor

  # Run on a Heme sample with custom sample ID
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_comp --sample-id MySample

  # Run on a Heme subpanel (MPN)
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel heme_mpn

  # Output to a specific directory
  python3 picard_qcsum.py --cram /path/to/sample.cram --panel pcp_normal --output-dir /path/to/output

  # List available panels
  python3 picard_qcsum.py --list-panels

  # Merge multiple qcsum files
  python3 picard_qcsum.py --merge sample1.qcsum.txt sample2.qcsum.txt -o /path/to/output

  # Batch process multiple CRAMs and merge results
  python3 picard_qcsum.py --batch sample1.cram sample2.cram --panel pcp_tumor -o /path/to/output
        """
    )

    parser.add_argument("--cram", "-c", help="Path to input CRAM/BAM file")
    parser.add_argument("--batch", "-b", nargs="+", metavar="CRAM", help="Batch process multiple CRAMs")
    parser.add_argument("--panel", "-p", choices=QCSumConfig.VALID_PANELS, help="Panel for QC cutoffs")
    parser.add_argument("--sample-id", "-s", help="Sample identifier (default: from filename)")
    parser.add_argument("--output-dir", "-o", help="Output directory")
    parser.add_argument("--reference", "-r", default=PICARD_REF, help="Reference FASTA")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Print commands only")
    parser.add_argument("--list-panels", action="store_true", help="List available panels")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--merge", "-m", nargs="+", metavar="FILE", help="Merge qcsum files")
    parser.add_argument("--merge-output", default="merged_qcsum.csv", help="Merged output filename")

    return parser.parse_args()
	

def main():
    """Main entry point."""
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.list_panels:
        list_panels()
        return 0

    # Handle merge mode
    if args.merge:
        try:
            output_file = merge_qcsum_files(
                qcsum_files=args.merge,
                output_dir=args.output_dir,
                output_name=args.merge_output
            )
            print(f"Merged qcsum written to: {output_file}")
            return 0
        except ValueError as e:
            logging.error(str(e))
            return 1

    # Handle batch mode - process multiple CRAMs and merge
    if args.batch:
        if not args.panel:
            print("Error: --panel is required for batch mode")
            return 1

        qcsum_files = []
        failed = []

        for cram_file in args.batch:
            try:
                result = run_qcsum(
                    cram_file=cram_file,
                    panel=args.panel,
                    output_dir=args.output_dir,
                    reference=args.reference,
                    dry_run=args.dry_run,
                )
                if not args.dry_run:
                    qcsum_files.append(result['qcsum_file'])
                    logging.info(f"Processed: {result['sample_id']} - Align:{result['alignment_qc']} Cov:{result['coverage_qc']}")
            except Exception as e:
                logging.error(f"Failed to process {cram_file}: {e}")
                failed.append(cram_file)

        if args.dry_run:
            return 0

        # Merge all qcsum files
        if qcsum_files:
            merged_file = merge_qcsum_files(
                qcsum_files=qcsum_files,
                output_dir=args.output_dir,
                output_name=args.merge_output
            )
            print(f"\nProcessed {len(qcsum_files)} samples, {len(failed)} failed")
            print(f"Merged qcsum: {merged_file}")

        return 1 if failed else 0

    # Validate required arguments for single-sample mode
    if not args.cram:
        print("Error: --cram or --batch is required (use --list-panels to see available panels)")
        return 1

    if not args.panel:
        print("Error: --panel is required")
        print("Use --list-panels to see available panels")
        return 1

    try:
        result = run_qcsum(
            cram_file=args.cram,
            panel=args.panel,
            sample_id=args.sample_id,
            output_dir=args.output_dir,
            reference=args.reference,
            dry_run=args.dry_run,
        )

        # Print summary
        if not args.dry_run:
            print("\n" + "=" * 60)
            print(f"Sample:       {result['sample_id']}")
            print(f"Panel:        {result['panel']}")
            print(f"Alignment QC: {result['alignment_qc']}")
            print(f"Coverage QC:  {result['coverage_qc']}")
            print(f"HSM file:     {result['hsm_file']}")
            print(f"QCSum file:   {result['qcsum_file']}")
            print("=" * 60)

        return 0

    except FileNotFoundError as e:
        logging.error(str(e))
        return 1
    except ValueError as e:
        logging.error(str(e))
        return 1
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        return e.returncode


if __name__ == "__main__":
    sys.exit(main())
