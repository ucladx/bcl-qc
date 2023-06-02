import sys
from helpers import *
from subprocess import call

BAM_OUTPUT_ROOT = "/mnt/pns/bams"
FASTQ_OUTPUT_ROOT = "/staging/hot/reads"
DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"

MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

class RunInfo:
    # More attributes can be added as needed, and will be accessible by all passes.
    def __init__(self):
        self.run_path = get_run_path()
        self.run_id = get_run_id(self.run_path)
        self.indices = get_indices(self.run_path)
        self.bed_path = DEFAULT_BED_PATH
        custom_passes = get_custom_passes()
        self.passes = custom_passes if custom_passes else MAIN_PASSES
        self.fastq_dir = f"{FASTQ_OUTPUT_ROOT}/{self.run_id}"

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
    for idx in run_info.indices:
        exec_pass("demux", run_info.run_path, idx, run_info.fastq_dir)

def align_pass(run_info):
    run_id = run_info.run_id
    for idx in run_info.indices:
        fastq_list = f"{run_info.fastq_dir}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            output_dir = f"{BAM_OUTPUT_ROOT}/{run_id}/{sample_id}"
            exec_pass("align", fastq_list, output_dir, run_info.bed_path, sample_id)

def multiqc_pass(run_info):
    save_occ_pf_plot(run_info.run_path)
    call(["bash", f"scripts/multiqc.sh", run_info.run_id])

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

def bclqc_run():
    run_info = RunInfo()
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    bclqc_run()
