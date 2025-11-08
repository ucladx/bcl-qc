import os
import re
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32
import shlex  # For safely constructing shell commands
import logging  # For more robust logging
from multiprocessing import Pool
import yaml

# --- Configuration ---
BAMS_DIR_DEFAULT = "/mnt/pns/bams/"
FASTQS_DIR_DEFAULT = "/staging/hot/reads/"

# Panel-specific config files
QCSUM_CONFIG_YAML = "config/qcsum_config.yaml"

PANEL_CONFIG = {
    "PCP": {
        "BED": "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed",
        "HUMAN_REF": "/staging/human/reference/hg38_alt_masked_graph_v3",
    },
    "HEME": {
        "BED": "/mnt/pns/tracks/goal_ucla_heme_221_exon_targets.hg38.bed",
        "HUMAN_REF": "/staging/human/reference/hg38_alt_masked_graph_v3",
    }
}

DEF_STEPS = [
    "demux",
    "align",
    "qc",
]

SAMPLEINFO_PANEL_TO_QCSUM_PANEL = {
    "Comprehensive Heme Panel": "heme_comp",
    "MPN Screen Panel (JAK2, CALR, MPL)": "heme_mpn",
    "JAK2": "heme_jak2",
    "Peripheral Blood Lymphoma Panel": "heme_pblp",
    "Tumor": "pcp_tumor",
    "Normal": "pcp_normal",
}

QC_SUM_HEADER = (
    "Sample,Sequencing_Platform,Pipeline_version,Alignment_QC,Coverage_QC,"
    "Total_Reads,%Reads_Aligned,Capture,Avg_Capture_Coverage,%On/Near_Bait_Bases,"
    "%On_Bait_Bases,FOLD_80_BASE_PENALTY,Avg_ROI_Coverage,MEDIAN_ROI_COVERAGE,"
    "MAX_ROI_COVERAGE,%ROI_1x,%ROI_20x,%ROI_100x,%ROI_250x,%ROI_500x"
)

PICARD_REF = "/mnt/pns/tracks/ref/hg38.fa"

class QCSumInfo:
    def __init__(self, panel, sample_id):
        self.panel = SAMPLEINFO_PANEL_TO_QCSUM_PANEL.get(panel, panel)
        self.sample_id = sample_id
        self.config = self.get_config()

    def get_config(self):
        with open(QCSUM_CONFIG_YAML) as f:
            yaml_obj = yaml.safe_load(f)
            if not yaml_obj:
                raise ValueError(f"Failed to load {QCSUM_CONFIG_YAML}. Ensure it is a valid YAML file.")
            if self.panel not in yaml_obj:
                raise ValueError(f"Panel '{self.panel}' not found in {QCSUM_CONFIG_YAML}.")
            if self.panel.startswith("heme_") and self.panel != "heme_comp":
                yaml_dict = yaml_obj.get("heme_comp", {})
                subpanel_dict = yaml_obj.get(self.panel, {})
                # override target intervals using subpanel BED
                target_intervals = subpanel_dict.get("target_intervals")
                if target_intervals:
                    yaml_dict["target_intervals"] = target_intervals
            else:
                yaml_dict = yaml_obj.get(self.panel, {})
            return {k: str(v) for k, v in yaml_dict.items()}

def parse_arguments():
    """
    Parses command line arguments for the bcl-qc pipeline.
    """
    parser = argparse.ArgumentParser(description="BCL QC pipeline for NovaSeq runs.")
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("--run-dir", required=True, help="Directory containing the run data to be processed")
    parser.add_argument("--fastqs-dir", help="Parent directory where FASTQs will be output", default=FASTQS_DIR_DEFAULT)
    parser.add_argument("--bams-dir", help="Parent directory where BAMs will be output", default=BAMS_DIR_DEFAULT)
    parser.add_argument("--sampleinfo", help="Path to the sampleinfo file for determining QC parameters.")
    parser.add_argument("--steps", nargs='+', help="Run only the specified steps (demux, align, qc)", default=DEF_STEPS)
    return parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def exec_command(cmd, executable=None, input_data=None):
    """
    Wrapper around subprocess.run, with enhanced error handling and logging.

    Args:
        cmd (list): Command to execute as a list of strings.
        executable (str, optional): Path to the executable if needed (e.g., 'sh' for shell scripts). Defaults to None.
        input_data (bytes, optional): Input data to pipe to the command's stdin. Defaults to None.

    Returns:
        subprocess.CompletedProcess: Result object from subprocess.run.

    Raises:
        subprocess.CalledProcessError: If the command fails (non-zero exit code).
    """
    cmd_str = ' '.join(map(shlex.quote, cmd)) # Quote each argument for safety
    logging.info(f"Running command: {cmd_str}")
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
            executable=executable,
            input=input_data # Handle input data for piping
        )
        if result.stderr:
            logging.warning(f"stderr: {result.stderr.decode()}") # Log stderr as warning, not always an error
        return result
    except subprocess.CalledProcessError as e:
        error_message = f"Command failed with return code: {e.returncode}\nCommand: {cmd_str}"
        if e.stderr:
            error_message += f"\nStderr: {e.stderr.decode()}" # Decode stderr for readability
        logging.error(error_message) # Log full error message
        raise # Re-raise the exception to stop pipeline execution

def demux(run_dir, samplesheet, fastq_output):
    """
    Perform demultiplexing on a NovaSeq run using Dragen BCLConvert.

    Args:
        run_dir (str): Path to the folder containing the NovaSeq run data.
        samplesheet (str): Path to the SampleSheet.csv to be used for this demultiplexing.
        fastq_output (str): Path to the folder where the demultiplexed FASTQs will be written.
    """
    logging.info(f"Starting demultiplexing for run directory: {run_dir}, samplesheet: {samplesheet}, output to: {fastq_output}")
    demux_cmd = [
        "dragen",
        "--bcl-conversion-only", "true",
        "--bcl-use-hw", "false",
        "--bcl-only-matched-reads", "true",
        "--bcl-input-directory", run_dir,
        "--sample-sheet", samplesheet,
        "--output-directory", fastq_output
    ]
    exec_command(demux_cmd)
    logging.info(f"Demultiplexing completed for output: {fastq_output}")

def align(fastq_list, bam_output, sample_id, panel="PCP"):
    """
    Perform alignment and variant calling using Dragen.

    Args:
        fastq_list (str): Path to the fastq list CSV file.
        bam_output (str): Path to the output directory for BAM/CRAM files.
        sample_id (str): Sample identifier.
        panel (str, optional): Sequencing panel type ('PCP' or 'HEME'). Defaults to "PCP".
    """
    logging.info(f"Starting alignment for sample: {sample_id}, fastq list: {fastq_list}, output to: {bam_output}, panel: {panel}")
    config = PANEL_CONFIG.get(panel)
    bed_file = config.get("BED")
    human_ref = config.get("HUMAN_REF")
    os.makedirs(bam_output, exist_ok=True)

    dragen_command = [
        "dragen",
        "--intermediate-results-dir", "/staging/tmp",
        "--enable-map-align", "true",
        "--enable-map-align-output", "true",
        "--output-format", "CRAM",
        "--generate-sa-tags", "true",
        "--enable-sort", "true",
        "--soft-read-trimmers", "polyg,quality",
        "--trim-min-quality", "2",
        "--ref-dir", human_ref,
        "--qc-coverage-tag-1", "target_bed",
        "--qc-coverage-region-1", bed_file,
        "--qc-coverage-reports-1", "cov_report",
        "--qc-coverage-ignore-overlaps", "true",
        "--enable-variant-caller", "true",
        "--vc-combine-phased-variants-distance", "6",
        "--vc-emit-ref-confidence", "GVCF",
        "--enable-hla", "true",
        "--fastq-list", fastq_list,
        "--fastq-list-sample-id", sample_id,
        "--output-directory", bam_output,
        "--output-file-prefix", sample_id,
    ]

    if panel == "HEME":
        dragen_command.extend([
            "--umi-enable", "true",
            "--umi-source", "qname",
            "--umi-correction-scheme", "random",
            "--umi-min-supporting-reads", "1",
            "--umi-metrics-interval-file", bed_file,
            "--vc-enable-umi-germline", "true",
            "--vc-enable-high-sensitivity-mode", "true",
        ])
    else:
        dragen_command.extend([
            "--enable-duplicate-marking", "true",
        ])

    exec_command(dragen_command)
    logging.info(f"Alignment completed for sample: {sample_id}, output: {bam_output}")

def multiqc(fastqs_dir, bams_dir):
    """
    Run MultiQC to aggregate QC reports.

    Args:
        fastqs_dir (str): Directory containing FASTQ files.
        bams_dir (str): Directory containing BAM/CRAM files.
    """
    logging.info(f"Starting MultiQC report generation, FASTQ dir: {fastqs_dir}, BAM dir: {bams_dir}")
    rm_cmd = ["rm", "-f", f"{bams_dir}/*/*.wgs_*.csv"] # Remove wgs coverage reports since they are irrelevant to targeted panels
    exec_command(rm_cmd)

    multiqc_cmd = [
        "multiqc",
        "--force",
        "--config", "config/multiqc_config.yaml", # Assuming config file exists at this path
        "--outdir", bams_dir,
        bams_dir,
        fastqs_dir,
    ]
    exec_command(multiqc_cmd)
    logging.info(f"MultiQC report generation completed in: {bams_dir}")

def qcsum_command(bam_file, sample, sample_dir, panel):
    """
    Run Picard CollectHsMetrics for QC summary.

    Args:
        bam_file (str): Path to the BAM/CRAM file.
        sample (str): Sample identifier.
        sample_dir (str): Directory for sample-specific outputs.
        panel (str): Sequencing panel type, based on sampleinfo.
    """
    logging.info(f"Running Picard CollectHsMetrics for sample: {sample}, panel: {panel}")
    qcsum_info = QCSumInfo(panel, sample).config

    bait_intervals = qcsum_info.get("bait_intervals")
    target_intervals = qcsum_info.get("target_intervals")
    output_file = os.path.join(sample_dir, f"{sample}.hsm.txt")

    picard_cmd = [
        "picard", "CollectHsMetrics",
        "I=" + bam_file,
        "O=" + output_file,
        "R=" + PICARD_REF,
        "BAIT_INTERVALS=" + bait_intervals,
        "TARGET_INTERVALS=" + target_intervals
    ]

    if os.path.exists(output_file):
        logging.info(f"Skipping Picard CollectHsMetrics for {sample} as output file already exists: {output_file}")
    else:
        exec_command(picard_cmd)
        logging.info(f"Picard CollectHsMetrics completed for sample: {sample}, output: {output_file}")

    perl_cmd = [
        "perl", "qcsum_metrics.pl",
        "--prefix", sample,
        "--qcfolder", sample_dir,
        "--pipeline_version", qcsum_info.get("pipeline_version"),
        "--platform", qcsum_info.get("platform"),
        "--pass_min_align_pct", qcsum_info.get("pass_min_align_pct"),
        "--fail_min_align_pct", qcsum_info.get("fail_min_align_pct"),
        "--covered", qcsum_info.get("covered"),
        "--pass_min_roi_pct", qcsum_info.get("pass_min_roi_pct"),
        "--fail_min_roi_pct", qcsum_info.get("fail_min_roi_pct"),
        "--pass_min_avgcov", qcsum_info.get("pass_min_avgcov"),
        "--fail_min_avgcov", qcsum_info.get("fail_min_avgcov"),
        "--pass_min_reads", qcsum_info.get("pass_min_reads"),
        "--fail_min_reads", qcsum_info.get("fail_min_reads"),
        "--capture", qcsum_info.get("capture"),
        "--capture_version", qcsum_info.get("capture_version")
    ]

    exec_command(perl_cmd)
    logging.info(f"qcsum_metrics.pl script executed for sample: {sample}, output: {output_file}")

def qcsum(bams_dir, sampleinfo):
    """
    Run qcsum.sh script for QC summary.

    Args:
        bams_dir (str): Directory containing BAM/CRAM files.
    """
    logging.info(f"Starting qcsum for BAM directory: {bams_dir}")
    sampleinfo_df = pd.read_csv(sampleinfo, sep="\t")
    qcsum_cmd_args = []
    qcsum_files = []
    for _, row in sampleinfo_df.iterrows():
        sample_id = row['Samples']
        panel = row.get('Panel', row.get('Tumor')) # "Panel" for Heme, "Tumor" for PCP
        bam_path = row.get('BAM Path')  # Use 'BAM Path' column from sampleinfo
        output_dir = os.path.join(bams_dir, sample_id)
        if not os.path.exists(output_dir): # output qc metrics into this run's directory even if sample is from a different run
            os.makedirs(output_dir, exist_ok=True)
        qcsum_cmd_args.append((bam_path, sample_id, output_dir, panel))
        qcsum_files.append(os.path.join(output_dir, f"{sample_id}.qcsum.txt"))
    with Pool() as pool:
        pool.starmap(qcsum_command, qcsum_cmd_args)
    with open(f"{bams_dir}/qcsum_mqc.csv", "w") as outfile:
        outfile.write(QC_SUM_HEADER + "\n")
        for qcsum_file in qcsum_files:
            with open(qcsum_file) as infile:
                lines = infile.readlines()
                outfile.writelines(lines[1:])  # Skip the first line (header) and write the QC info
        outfile.write("\n")
    logging.info(f"qcsum.sh execution completed for directory: {bams_dir}")

def demux_samples(run_dir, fastqs_dir):
    """
    Execute the demultiplexing step.

    Args:
        run_dir (str): Path to the run directory.
        fastqs_dir (str): Path to the FASTQ output directory.
    """
    logging.info("Starting demux step")
    if os.path.exists(fastqs_dir):
        error_msg = "Skipping demux as output directory already exists: " + fastqs_dir
        logging.error(error_msg)
        raise Exception(error_msg)
    os.makedirs(fastqs_dir, exist_ok=True)
    samplesheets = get_samplesheets(run_dir)
    for samplesheet in samplesheets:
        fastq_output = f"{fastqs_dir}/{get_index(samplesheet)}"
        demux(run_dir, samplesheet, fastq_output)
    logging.info("Demux step completed")

def align_samples(fastqs_dir, bams_dir):
    """
    Execute the alignment step.

    Args:
        fastqs_dir (str): Path to the FASTQ directory.
        bams_dir (str): Path to the BAM output directory.
    """
    logging.info("Starting align step")

    fastq_lists = [os.path.join(root, file) for root, _, files in os.walk(fastqs_dir) for file in files if file == "fastq_list.csv"]
    for fastq_list in fastq_lists:
        if "/U5N2_I10/" in fastq_list:
            panel = "HEME"
        elif "/I10/" in fastq_list:
            panel = "PCP"
        elif "/I8N2_N10/" in fastq_list:
            logging.info(f"Skipping alignment for exome samples in: {fastq_list}") # Skip exome samples
            continue
        else:
            error_msg = "Could not determine bed file for fastq list: " + fastq_list
            logging.error(error_msg)
            raise Exception(error_msg)
        for sample_id in get_sample_ids(fastq_list):
            bam_output = os.path.join(bams_dir, sample_id)
            if os.path.exists(bam_output):
                logging.info(f"Skipping alignment for {sample_id} as output directory already exists: {bam_output}")
                continue
            align(fastq_list, bam_output, sample_id, panel)
    logging.info("Align step completed")

def qc_samples(run_dir, fastqs_dir, bams_dir, sampleinfo):
    """
    Execute the QC step, including saving plots, running MultiQC and qcsum.

    Args:
        run_dir (str): Path to the run directory.
        fastqs_dir (str): Path to the FASTQ directory.
        bams_dir (str): Path to the BAM directory.
        sampleinfo (str): Path to sampleinfo with which to determine which BED/ilists to use for qcsum
    """
    logging.info("Starting qc step")
    save_occ_pf_plot(run_dir, bams_dir)
    qcsum(bams_dir, sampleinfo)
    multiqc(fastqs_dir, bams_dir)
    logging.info("QC step completed")

def get_index(file_name):
    """
    Extract the index from a samplesheet filename.

    Args:
        file_name (str): Samplesheet filename (e.g., "SampleSheet_I10.csv").

    Returns:
        str: The index string (e.g., "I10").
    """
    # Assumes samplesheet format is "SampleSheet_{index}.csv"
    return file_name.split('/')[-1].replace("SampleSheet_", "").replace(".csv", "")

def get_samplesheets(run_dir):
    """
    Get a list of samplesheet filenames in a directory.

    Args:
        run_dir (str): Path to the run directory.

    Returns:
        list: List of samplesheet filenames.

    Raises:
        Exception: If no samplesheets are found in the directory.
    """
    samplesheets = [os.path.join(run_dir, f) for f in os.listdir(run_dir) if "SampleSheet_" in f]
    if not samplesheets:
        error_msg = "Error: No samplesheets found in " + run_dir
        logging.error(error_msg)
        raise Exception(error_msg)
    return samplesheets

def get_sample_ids(fastq_list):
    """
    Extract sample IDs from a fastq list CSV file.

    Args:
        fastq_list (str): Path to the fastq_list.csv file.

    Returns:
        set: Set of unique sample IDs (RGSM column).
    """
    fastq_list_df = pd.read_csv(fastq_list)
    return set(fastq_list_df['RGSM'])

def parse_run_metrics(run_dir):
    """
    Parses Illumina Interop run metrics from a run directory.

    Args:
        run_dir (str): Path to the Illumina run directory.

    Returns:
        pd.DataFrame or None: DataFrame of interop imaging table if data is found, None otherwise.
    """
    logging.info(f"Parsing run metrics from: {run_dir}")
    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    try:
        run_metrics.read(run_dir, valid_to_load)
    except Exception as e: # Catch broad exception for interop errors
        logging.error(f"Error occurred trying to open {run_dir}: {e}")
        return None

    # Set up data table
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    ncolumns = columns.size()
    if ncolumns == 0:
        logging.warning("No interop data found in run metrics.")
        return None # no data

    headers = []
    for i in range(ncolumns):
        column = columns[i]
        if column.has_children():
            headers.extend(
                [f"{column.name()} ({subname})" for subname in column.subcolumns()])
        else:
            headers.append(column.name())

    column_count = py_interop_table.count_table_columns(columns)
    row_offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, row_offsets)
    data = zeros((row_offsets.size(), column_count), dtype=float32)

    # Populate table
    py_interop_table.populate_imaging_table_data(
        run_metrics, columns, row_offsets, data.ravel()
    )

    logging.info(f"Interop run metrics parsed successfully from: {run_dir}")
    return pd.DataFrame(data, columns=headers)

def save_occ_pf_plot(run_dir, output_dir):
    """
    Saves a % Occupied x % Pass Filter scatter plot to the output directory.

    Args:
        run_dir (str): Path to the Illumina run directory.
        output_dir (str): Path to the output directory to save the plot.
    """
    logging.info(f"Generating and saving Occupied vs Pass Filter plot for: {run_dir}, output to: {output_dir}")
    df = parse_run_metrics(run_dir)
    if df is None:
        logging.warning("Unable to parse Interop files, skipping plot generation.")
        return

    x = "% Occupied"
    y = "% Pass Filter"
    views = ["Lane"] # can add "Tile" or "Cycle" for more views
    for view in views:
        sns.scatterplot(
            data=df,
            x=x,
            y=y,
            hue=view,
            alpha=0.5,
            s=8,
        )
        plt.xlim([0, 100])
        plt.ylim([50, 100])
        plt.legend(title=view, bbox_to_anchor=[1.2, 0.9])
        plt.tight_layout()

        image_path = os.path.join(output_dir, f"occ_pf_{view.lower()}_mqc.jpg")
        logging.info(f"Saving Occupied vs Pass Filter graph to: {image_path}")
        plt.savefig(image_path, dpi=300)
        plt.close()
    logging.info(f"Occupied vs Pass Filter plot saved to: {output_dir}")

def bclqc_run():
    """
    Main function to execute the bcl-qc pipeline.
    """
    args = parse_arguments()
    steps = args.steps
    run_dir = args.run_dir
    sampleinfo = args.sampleinfo
    run_name = run_dir.split('/')[4]
    fastqs_dir = args.fastqs_dir + '/' + run_name
    bams_dir = args.bams_dir + '/' + run_name

    if "demux" in steps:
        demux_samples(run_dir, fastqs_dir)
    if "align" in steps:
        align_samples(fastqs_dir, bams_dir)
    if "qc" in steps:
        if not sampleinfo:
            logging.error("Sampleinfo file not found, cannot determine panel for QCSum.")
        qc_samples(run_dir, fastqs_dir, bams_dir, sampleinfo)

    logging.info("BCL QC pipeline run completed.")

if __name__ == "__main__":
    bclqc_run()
