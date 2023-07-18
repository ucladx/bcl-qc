import sys
from helpers import *
from subprocess import call

BAM_OUTPUT_ROOT = "/mnt/pns/bams"
FASTQ_OUTPUT_ROOT = "/staging/hot/reads"
DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"
RUN_DIR = "/mnt/pns/runs/"
HUMAN_REF = "/staging/human/reference/hg38_alt_masked_graph_v2"

MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

def azure_pass(run_info):
    print("Azure upload placeholder...")
    return

def megaqc_pass(run_info):
    print("MegaQC placeholder...")
    return

def demux_cmd(run_path, samplesheet, fastq_output):
    call(["dragen",
          "--bcl-conversion-only", "true",
          "--bcl-use-hw", "false",
          "--bcl-only-matched-reads", "true",
          "--bcl-input-directory", run_path,
          "--sample-sheet", samplesheet,
          "--output-directory", fastq_output])

def demux_pass(run_info):
    run_name = run_info.run_name
    call(["mkdir", "-p", f"/staging/hot/reads/{run_name}/"])
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_name}/SampleSheet_{idx}.csv"
        fastq_output = f"/staging/hot/reads/{run_name}/{idx}"
        demux_cmd(run_info.run_path, samplesheet, fastq_output)

def align_cmd(fastq_list, bam_output, bed_path, sample_id):
    call(["mkdir", "-p", bam_output])
    call(["dragen",
          "--enable-map-align", "true",
          "--enable-map-align-output", "true",
          "--output-format", "BAM",
          "--enable-duplicate-marking", "true",
          "--generate-sa-tags", "true",
          "--enable-sort", "true",
          "--soft-read-trimmers", "polyg,quality",
          "--trim-min-quality", "2",
          "--ref-dir", HUMAN_REF,
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
          "--output-file-prefix", sample_id,
    ])

def align_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_name}/{sample_id}"
            align_cmd(fastq_list, bam_output, run_info.bed_path, sample_id)

def multiqc_cmd(run_name):
    call(["rm", "-f", "/mnt/pns/bams/$run_name/*/*.wgs_*.csv"])
    call(["multiqc",
          "--force",
          "--config", "config/multiqc_config.yaml",
          "--outdir", f"/mnt/pns/bams/{run_name}/",
           # dirs to scan
          f"/mnt/pns/bams/{run_name}",
          f"/staging/hot/reads/{run_name}/",
    ])


def multiqc_pass(run_info):
    save_occ_pf_plot(run_info.run_path, run_info.run_name)
    multiqc_cmd(run_info.run_name)

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
    def __init__(self, run_path, args):
        self.run_path = run_path
        self.run_name = get_run_name(run_path)
        self.bed_path = DEFAULT_BED_PATH
        self.indices = get_indices(self.run_path)
        custom_passes = args.get('P')
        self.passes = custom_passes if custom_passes else MAIN_PASSES

def bclqc_run(run_path, args=None):
    run_info = RunInfo(run_path, args)
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    run_path = get_run_path()
    args = get_args()
    bclqc_run(run_path, args)
