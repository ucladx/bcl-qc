import sys
from helpers import *
from subprocess import call

HELP_MSG = """BCLQC
Usage: python3 bclqc.py [options] fastqs_dir bams_dir run_dir

fastqs_dir is the parent directory where FASTQs will be output
bams_dir is the parent directory where BAMs will be output
run_dir is the directory containing the run data to be processed

For example:
`python3 bclqc.py [options] /fastqs /bams /runs/210930`

In this case /fastqs/210930 and /bams/210930 will be created if they don't exist,
and FASTQs and BAMs will be written to those directories respectively.

Options:
    -P [pass1,pass2,...]    Run only the specified passes
    -h, --help              Print this help message
    -d, --dry-run           Dry run: print commands without executing them"""

DEFAULT_BED_FILE = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"
DEFAULT_HUMAN_REF = "/staging/human/reference/hg38_alt_masked_graph_v2"

MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

def demux(run_dir, samplesheet, fastq_output,
          exec_cmd=call):
    """
    Perform demultiplexing on a NovaSeq run using Dragen BCLConvert
    Inputs:
        run_dir: path to the folder containing the NovaSeq run data
        samplesheet: path to the SampleSheet.csv to be used for this demultiplexing
        fastq_output: path to the folder where the demultiplexed FASTQs will be written
    """
    exec_cmd(["dragen",
          "--bcl-conversion-only", "true",
          "--bcl-use-hw", "false",
          "--bcl-only-matched-reads", "true",
          "--bcl-input-directory", run_dir,
          "--sample-sheet", samplesheet,
          "--output-directory", fastq_output])

def align(fastq_list, bam_output, bed_file, sample_id, exec_cmd=call):
    exec_cmd(["mkdir", "-p", bam_output])
    exec_cmd(["dragen",
          "--enable-map-align", "true",
          "--enable-map-align-output", "true",
          "--output-format", "BAM",
          "--enable-duplicate-marking", "true",
          "--generate-sa-tags", "true",
          "--enable-sort", "true",
          "--soft-read-trimmers", "polyg,quality",
          "--trim-min-quality", "2",
          "--ref-dir", DEFAULT_HUMAN_REF,
          "--intermediate-results-dir", "/staging/tmp",
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
    ])

def multiqc_cmd(fastqs_dir, bams_dir, exec_cmd=call):
    exec_cmd(["rm", "-f", f"{bams_dir}/*/*.wgs_*.csv"])
    exec_cmd(["multiqc",
          "--force",
          "--config", "config/multiqc_config.yaml",
          "--outdir", bams_dir,
           # dirs to scan
          bams_dir,
          fastqs_dir,
    ])

def demux_pass(run_info):
    run_dir = run_info.run_dir
    fastqs_dir = run_info.fastqs_dir
    exec_cmd = run_info.exec_cmd
    exec_cmd(["mkdir", "-p", fastqs_dir])
    for idx in run_info.indices:
        samplesheet = f"{run_dir}/SampleSheet_{idx}.csv"
        fastq_output = f"{fastqs_dir}/{idx}"
        demux(run_dir, samplesheet, fastq_output, exec_cmd)

def align_pass(run_info):
    fastqs_dir = run_info.fastqs_dir
    bams_dir = run_info.bams_dir
    for idx in run_info.indices:
        fastq_list = f"{fastqs_dir}/{idx}/Reports/fastq_list.csv"
        if not os.path.exists(fastq_list):
            raise Exception(f"Alignment Error: {fastq_list} does not exist")
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"{bams_dir}/{sample_id}"
            align(fastq_list, bam_output, run_info.bed_file, sample_id, run_info.exec_cmd)

def multiqc_pass(run_info):
    save_occ_pf_plot(run_info.run_dir, run_info.bams_dir, run_info.exec_cmd)
    multiqc_cmd(run_info.fastqs_dir, run_info.bams_dir, run_info.exec_cmd)

def get_pass_f(pass_name):
    pass_f_name = pass_name + "_pass"
    # look into global namespace for pass function
    for name, f in globals().items():
        if name == pass_f_name:
            if callable(f) and f.__module__ == __name__:
                return f
    return None

def execute_pass(pass_name, run_info):
    pass_function = get_pass_f(pass_name)
    if pass_function: # if a custom pass function is defined, call it first
        print(f"Running {pass_name}...")
        pass_function(run_info)
    else:
        print(f"Pass not found: {pass_name}\n"
              f"Define a function called {pass_name}_pass in bcl-qc.py")

class RunInfo:
    def __init__(self):
        [fastqs_dir, bams_dir, run_dir] = get_arg_dirs()
        self.run_dir = run_dir
        self.bams_dir = bams_dir
        self.fastqs_dir = fastqs_dir
        # TODO allow BED paths to be specified per sample in samplesheet
        self.bed_file = DEFAULT_BED_FILE
        self.indices = get_indices(self.run_dir)
        self.passes = get_passes()
        self.exec_cmd = get_exec_cmd()

def bclqc_run():
    if "-h" in sys.argv or "--help" in sys.argv:
        print(HELP_MSG)
        return
    run_info = RunInfo()
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    bclqc_run()
