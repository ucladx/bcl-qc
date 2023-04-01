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

DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"

# Main passes are run by default
MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

def azure_pass(run_info):
    print("TEST Azure Upload...")
    return

def megaqc_pass(run_info):
    print("TEST MegaQC run...")

def demux_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_name}/SampleSheet_{idx}.csv"
        fastq_output = f"/staging/hot/reads/{run_name}/"
        call(["bash", "scripts/demux.sh",
                run_info.run_path,
                samplesheet,
                fastq_output])

def align_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_name}/$sample_id"
            call(["bash", "scripts/align.sh",
                    fastq_list,
                    bam_output,
                    run_info.bed_path,
                    sample_id])

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
    save_occ_pf_plot(run_info.run_path)
    call(["bash", f"scripts/multiqc.sh", run_info.run_name])

def get_pass_f(pass_name):
    pass_f_name = pass_name + "_pass"
    for name, f in globals().items():
        if name == pass_f_name:
            if callable(f) and f.__module__ == __name__:
                return f
    return None

def compute_passes(args):
    passes = MAIN_PASSES
    if args.get('a'): passes += "azure"
    if args.get('M'): passes += "megaqc"
    print("passes to run: ", passes)
    return passes

def execute_pass(pass_name, run_info):
    pass_function = get_pass_f(pass_name)
    if pass_function: # if a custom pass function is defined, call it first
        pass_function(run_info)
    else:
        print(f"I couldn't execute pass: {pass_name} because a pass function \n"
              f"Define a function called {pass_name}_pass in bcl-qc.py")

class RunInfo:
    def __init__(self, run_path, args):
        self.run_path = run_path
        self.indices = get_indices(run_path)
        self.run_name = get_run_name(run_path)
        self.bed_path = get_bed_path(args, DEFAULT_BED_PATH)
        custom_passes = args.get('O')
        self.passes = custom_passes if custom_passes else compute_passes(args)

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