import sys
from helpers import *
from subprocess import call

#   ____   _____ _           ____   _____ 
#  |  _ \ / ____| |         / __ \ / ____|
#  | |_) | |    | |  ______| |  | | |     
#  |  _ <| |    | | |______| |  | | |     
#  | |_) | |____| |____    | |__| | |____ 
#  |____/ \_____|______|    \___\_\\_____|
#


### USAGE ###
# Basic Usage:
# python bclqc.py <run_path>
#   this will run all MAIN_PASSES
# 
# Flags:
#   -P <pass_name> <pass_name> ...
#   run only the specified passes
#

DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"

MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

def shell_exec(*args):
    call(["bash", *args])

def exec_pass(script_name, *args):
    get_script = lambda script_name: f"scripts/{script_name}.sh"
    shell_exec(get_script(script_name), *args)

def azure_pass(run_info):
    print("TEST Azure Upload...")
    return

def megaqc_pass(run_info):
    print("TEST MegaQC run...")

def demux_pass(run_info):
    run_id = run_info.run_id
    call(["bash", "mkdir", f"/staging/hot/reads/{run_id}/"])
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_id}/SampleSheet_{idx}.csv"
        fastq_output = f"/staging/hot/reads/{run_id}/{idx}"
        exec_pass("demux", run_info.run_path, samplesheet, fastq_output)

def align_pass(run_info):
    run_id = run_info.run_id
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_id}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_id}/{sample_id}"
            call(["bash", "mkdir", bam_output])
            exec_pass("align", fastq_list, bam_output, run_info.bed_path, sample_id)

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
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
        pass_function(run_info)
    else:
        print(f"I couldn't execute pass: {pass_name} because a pass function \n"
              f"Define a function called {pass_name}_pass in bcl-qc.py")

class RunInfo:
    def __init__(self, run_path, args):
        self.run_path = run_path
        self.bed_path = DEFAULT_BED_PATH
        custom_passes = args.get('P')
        self.passes = custom_passes if custom_passes else MAIN_PASSES

    def get_run_id(self):
        return self.run_path.split('/')[-2]

def bclqc_run():
    run_path = sys.argv[-1]
    if run_path[-1] != '/':  run_path += '/'
    args = parse_args(sys.argv[1:-1])
    run_info = RunInfo(run_path, args)
    for pass_name in run_info.passes:
        execute_pass(pass_name, run_info)

if __name__ == "__main__":
    bclqc_run()