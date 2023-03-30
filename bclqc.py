import sys
from subprocess import call
from os.path import exists
from helpers import *

### TODO ###
# Better file organization
# Sample specific bed files (parse samplesheet)
# Support passes config file? 

# Unit tests:
#     sample sheet generation
#     run detection
#     Pass manipulation

DEFAULT_BED = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"
REF_PATH = "/staging/human/reference/hg38_alt_masked_graph_v2"

# Main passes are run by default
MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

class RunInfo:
    def __init__(self, run_path, args):
        self.run_path = run_path
        self.indices = get_indices(run_path)
        self.run_name = get_run_name(run_path)
        self.bed_path = get_bed_path(args, DEFAULT_BED)
        custom_passes = args.get('O')
        self.passes = custom_passes if custom_passes else compute_passes(args)

def demux_cmd(samplesheet, run_path, fastq_output):
    call([
        "dragen",
        "--bcl-conversion-only", "true",
        "--bcl-use-hw", "false",
        "--bcl-only-matched-reads", "true",
        "--bcl-input-directory", run_path,
        "--sample-sheet", samplesheet,
        "--output-directory", fastq_output,
    ])

def align_cmd(fastq_list, sample_id, bed_path, ref_path, bam_output):
    call([
        "dragen",
        "--enable-map-align", "true",
        "--enable-map-align-output", "true",
        "--output-format", "BAM",
        "--enable-duplicate-marking", "true",
        "--generate-sa-tags", "true",
        "--enable-sort", "true",
        "--soft-read-trimmers", "polyg,quality",
        "--trim-min-quality", "2",
        "--ref-dir", ref_path,
        "--intermediate-results-dir", "/staging/tmp",
        "--qc-coverage-tag-1", "target_bed",
        "--qc-coverage-region-1", bed_path,
        "--qc-coverage-reports-1", "cov_report",
        "--qc-coverage-ignore-overlaps", "true",
        "--enable-variant-caller", "true",
        "--vc-combine-phased-variants-distance", "6",
        "--vc-emit-ref-confidence", "GVCF",
        "--enable-hla", "true",
        "--fastq-list", fastq_list,
        "--fastq-list-sample-id", sample_id,
        "--output-directory", bam_output,
        "--output-file-prefix", sample_id
    ])

def multiqc_cmd(reads_input, run_path, bam_dir):
    save_occ_pf_plot(run_path)
    # Remove *.wgs_*.csv files
    for root, _, files in os.walk(bam_dir):
        for file in files:
            if ".wgs_" in file and file.endswith(".csv"):
                os.remove(os.path.join(root, file))
    call([
        "multiqc",
        "--force",
        "--config", "./config/multiqc_config.yaml",
        "--outdir", bam_dir,
        bam_dir,
        reads_input
    ])

def demux_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_name}/SampleSheet_{idx}.csv"
        fastq_output = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        demux_cmd(samplesheet, run_info.run_path, fastq_output)

def align_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_name}/{sample_id}"
            align_cmd(fastq_list, sample_id, run_info.bed_path, REF_PATH, bam_output)

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
    run_name = run_info.run_name
    bam_dir = f"/mnt/pns/bams/{run_name}"
    reads_input = f"/staging/hot/reads/{run_name}"
    multiqc_cmd(reads_input, run_info.run_path, bam_dir)

def azure_pass(run_info):
    print("TEST Azure Upload...")
    return

def megaqc_pass(run_info):
    print("TEST MegaQC run...")

def get_pass_function(pass_name):
    function_name = f"{pass_name}_pass"
    return globals().get(function_name)

def compute_passes(args):
    passes = MAIN_PASSES
    if args.get('a'): passes += "azure"
    if args.get('M'): passes += "megaqc"
    print("passes to run: ", passes)
    return passes

def execute_pass(pass_name, run_info):
    pass_function = get_pass_function(pass_name)
    if pass_function: # if a custom pass function is defined, call it first
        pass_function(run_info)
    else:
        print(f"I couldn't execute pass: {pass_name} because a pass function \n"
              f"Define a function called {pass_name}_pass in bcl-qc.py")

def bclqc_run(run_path, args):
    run_info = RunInfo(run_path, args)
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    # absolute path to a specific run's output directory
    run_path = sys.argv[-1]
    if run_path[-1] != '/':  run_path += '/'
    args = parse_args(sys.argv[1:-1])
    bclqc_run(run_path, args)