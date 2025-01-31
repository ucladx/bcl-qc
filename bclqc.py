import os
import re
from subprocess import call
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32

HUMAN_REF = "/staging/human/reference/hg38_alt_masked_graph_v3"

# Maps barcode indices to the bed file that corresponds to the target regions for that panel
# I10 = Pan-Cancer Panel
# U5N2_I10 = Heme Panel
BEDS = {
    "PCP": "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed",
    "HEME": "/mnt/pns/tracks/goal_ucla_heme_221_exon_targets.hg38.bed",
}

DEF_PASSES = [
    "demux",
    "align",
    "qc",
]

def demux(run_dir, samplesheet, fastq_output):
    """
    Perform demultiplexing on a NovaSeq run using Dragen BCLConvert
    Inputs:
        run_dir: path to the folder containing the NovaSeq run data
        samplesheet: path to the SampleSheet.csv to be used for this demultiplexing
        fastq_output: path to the folder where the demultiplexed FASTQs will be written
    """
    call(["dragen",
          "--bcl-conversion-only", "true",
          "--bcl-use-hw", "false",
          "--bcl-only-matched-reads", "true",
          "--bcl-input-directory", run_dir,
          "--sample-sheet", samplesheet,
          "--output-directory", fastq_output])

def align(fastq_list, bam_output, sample_id, panel="PCP"):
    bed_file = BEDS[panel]
    call(["mkdir", "-p", bam_output])
    dragen_command = [
        "dragen",
        "--intermediate-results-dir", "/staging/tmp",
        "--enable-map-align", "true",
        "--enable-map-align-output", "true",
        "--output-format", "BAM" if panel == "PCP" else "CRAM",
        "--generate-sa-tags", "true",
        "--enable-sort", "true",
        "--soft-read-trimmers", "polyg,quality",
        "--trim-min-quality", "2",
        "--ref-dir", HUMAN_REF,
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
        "--vc-output-evidence-bam", "true",
        "--vc-evidence-bam-output-haplotypes", "true",
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

    call(dragen_command)

def multiqc(fastqs_dir, bams_dir):
    call(["rm", "-f", f"{bams_dir}/*/*.wgs_*.csv"]) # Remove wgs coverage reports since we only care about panel coverage
    call(["multiqc",
          "--force",
          "--config", "config/multiqc_config.yaml",
          "--outdir", bams_dir,
           # dirs to scan
          bams_dir,
          fastqs_dir,
    ])

def qcsum(bams_dir):
    for sample in os.listdir(bams_dir):
        if not os.path.isfile(f"{bams_dir}/{sample}/{sample}.cram"):
            continue
        call(["qcsum.sh",
              f"{bams_dir}/{sample}/{sample}.bam",
              sample])

def demux_pass(run_dir, fastqs_dir):
    call(["mkdir", "-p", fastqs_dir])
    samplesheets = get_samplesheets(run_dir)
    for samplesheet in samplesheets:
        fastq_output = f"{fastqs_dir}/{get_index(samplesheet)}"
        demux(run_dir, samplesheet, fastq_output)

def align_pass(fastqs_dir, bams_dir):
    fastq_lists = [os.path.join(root, file) for root, _, files in os.walk(fastqs_dir) for file in files if file == "fastq_list.csv"]
    for fastq_list in fastq_lists:
        if "/U5N2_I10/" in fastq_list:
            panel = "HEME"
        elif "/I10/" in fastq_list:
            panel = "PCP"
        elif "/I8N2_N10/" in fastq_list:
            continue # skip aligning exome samples
        else:
            raise Exception("Could not determine bed file for fastq list: " + fastq_list)
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"{bams_dir}/{sample_id}"
            align(fastq_list, bam_output, sample_id, panel)

def qc_pass(run_dir, fastqs_dir, bams_dir):
    save_occ_pf_plot(run_dir, bams_dir)
    multiqc(fastqs_dir, bams_dir)
    qcsum(bams_dir)

def get_run_name(run_dir: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    split_fields = run_dir.split('/')
    if split_fields[-1] == '':
        return split_fields[-2]
    else:
        return split_fields[-1]

def is_samplesheet(file_name):
    return re.search("^SampleSheet_.*\.csv$", file_name)

def get_index(file_name):
    # assumes samplesheet format is "SampleSheet_{index}.csv"
    return file_name.replace("SampleSheet_", "").replace(".csv", "")

def get_samplesheets(run_dir):
    samplesheets = [f for f in os.listdir(run_dir) if is_samplesheet(f)]
    if not samplesheets:
        raise("Error: No samplesheets found in " + run_dir)
    return samplesheets

def get_sample_ids(fastq_list):
    fastq_list_df = pd.read_csv(fastq_list)
    return set(fastq_list_df['RGSM'])

def parse_run_metrics(run_dir: str):
    """
    Returns a DataFrame of `interop_imaging_table` or returns None if no data is found.

    More info on `interop_imaging_table` and the `interop` library:
    http://illumina.github.io/interop/index.html
    https://github.com/Illumina/interop
    """

    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    try:
        run_metrics.read(run_dir, valid_to_load)
    except:
        print(f"Error occured trying to open {run_dir}")
        return None

    # Set up data table
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    ncolumns = columns.size()
    if ncolumns == 0: return None # no data

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

    return pd.DataFrame(data, columns=headers)

def save_occ_pf_plot(run_dir, output_dir):
    """
    Saves a % Occupied x % Pass Filter scatter to `SAVE_DIR`.
    """
    df = parse_run_metrics(run_dir)
    if df is None:
        print("Unable to parse Interop files.")
        return
    x = "% Occupied"
    y = "% Pass Filter"
    views = ["Lane"] # can add "Tile" or "Cycle" for more views
    for view in views:
        scatterplot(
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

        image_path = output_dir + f"/occ_pf_{view.lower()}_mqc.jpg"
        print("saving occ pf graph to " + image_path)
        plt.savefig(image_path, dpi=300)
        plt.close()

parser = argparse.ArgumentParser(description="bcl-qc pipeline")
parser.add_argument("--run_dir", required=True, help="Directory containing the run data to be processed")
parser.add_argument("--fastqs_dir", required=True, help="Parent directory where FASTQs will be output")
parser.add_argument("--bams_dir", required=True, help="Parent directory where BAMs will be output")
parser.add_argument("--passes", nargs='+', help="Run only the specified passes", default=DEF_PASSES)
args = parser.parse_args()

def get_arg_dirs():
    dirs = [args.run_dir, args.fastqs_dir, args.bams_dir]
    for d in dirs:
        if not os.path.exists(d):
            raise Exception(f"Given directory {d} does not exist")
    return dirs

def bclqc_run():
    passes = args.passes
    run_dir, fastqs_dir, bams_dir = get_arg_dirs()
    if "demux" in passes:
        demux_pass(run_dir, fastqs_dir)
    if "align" in passes:
        align_pass(fastqs_dir, bams_dir)
    if "qc" in passes:
        qc_pass(fastqs_dir, bams_dir)

if __name__ == "__main__":
    bclqc_run()
