import sys
from helpers import *

BAM_OUTPUT_ROOT = "/mnt/pns/bams"
FASTQ_OUTPUT_ROOT = "/staging/hot/reads"
DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"
RUN_DIR = "/mnt/pns/runs/"

MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

class RunInfo:
    # More attributes can be added as needed, and will be accessible by all passes.
    def __init__(self, run_name, args):
        self.run_name = run_name
        self.run_path = RUN_DIR + run_name + "/"
        self.bed_path = DEFAULT_BED_PATH
        self.indices = get_indices(self.run_path)
        custom_passes = get_custom_passes()
        self.passes = custom_passes if custom_passes else MAIN_PASSES

    def display(self):
        print(f"Run Path: {self.run_path}")
        print(f"Run ID: {self.run_id}")
        print(f"Indices: {self.indices}")
        print(f"Bed Path: {self.bed_path}")
        print(f"Passes: {self.passes}")

def azure_pass(run_info):
    print("Azure upload placeholder...")
    return

def megaqc_pass(run_info):
    print("MegaQC placeholder...")
    return

def demux_pass(run_info):
    run_name = run_info.run_name
    shell_exec(["mkdir", f"/staging/hot/reads/{run_name}/"])
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_name}/SampleSheet_{idx}.csv"
        fastq_output = f"/staging/hot/reads/{run_name}/{idx}"
        exec_pass("demux", run_info.run_path, samplesheet, fastq_output)

def align_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_name}/{sample_id}"
            shell_exec(["mkdir", bam_output])
            exec_pass("align", fastq_list, bam_output, run_info.bed_path, sample_id)

def multiqc_pass(run_info):
    save_occ_pf_plot(run_info.run_path)
    exec_pass("multiqc", run_info.run_name)

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
        print(f"I couldn't execute pass: {pass_name} because a pass function \n"
              f"Define a function called {pass_name}_pass in bcl-qc.py")

class RunInfo:
    def __init__(self, run_name, args):
        self.run_name = run_name
        self.run_path = RUN_DIR + run_name + "/"
        self.bed_path = DEFAULT_BED_PATH
        self.indices = get_indices(self.run_path)
        custom_passes = args.get('P')
        self.passes = custom_passes if custom_passes else MAIN_PASSES

def bclqc_run(run_name, args=None):
    run_info = RunInfo(run_name, args)
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    run_name = sys.argv[-1]
    args = parse_args(sys.argv[1:-1])
    bclqc_run(run_name, args)
