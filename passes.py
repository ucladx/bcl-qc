from subprocess import call
from os.path import exists
from helpers import *

DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"

# Main passes are run by default
MAIN_PASSES = [
    "demux",
    "align",
    "multiqc",
]

# Aux passes are defined with a flag that can be given to enable them
AUX_PASSES = {
    'M' : "megaqc",
    'a' : "azure",
}

class RunInfo:
    def __init__(self, run_path, args):
        user_bed_path = get_flag_args('b', args)
        self.run_path = run_path
        self.indices = get_indices(run_path)
        self.run_name = get_run_name(run_path)
        self.bed_path = user_bed_path[0] if user_bed_path else DEFAULT_BED_PATH

def azure_pass(run_info):
    print("TEST Azure Upload\n-------------------------")
    return

def demux_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        samplesheet = f"/mnt/pns/runs/{run_name}/SampleSheet_{idx}.csv"
        fastq_dir = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        call(["bash", f"scripts/demux.sh",
                run_info.run_path,
                samplesheet,
                fastq_dir])

def align_pass(run_info):
    run_name = run_info.run_name
    for idx in run_info.indices:
        fastq_list = f"/staging/hot/reads/{run_name}/{idx}/Reports/fastq_list.csv"
        for sample_id in get_sample_ids(fastq_list):
            bam_output = f"/mnt/pns/bams/{run_name}/$sample_id"
            call(["bash", f"scripts/align.sh",
                    fastq_list,
                    bam_output,
                    run_info.bed_path,
                    sample_id])

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
    save_occ_pf_plot(run_info.run_path)
    call(["bash", f"scripts/multiqc.sh", run_info.run_name])

def handle_skip_flag(passes, skipped_passes):
    new_passes = passes.copy()
    for pass_name in skipped_passes:
        if pass_name.upper() == "ALL":
            return []
        elif pass_name.upper() == "MAIN":
            for pass_name in MAIN_PASSES: new_passes.remove(pass_name)
        else:
            new_passes.remove(pass_name)
    return new_passes

def get_custom_pass(pass_name):
    custom_pass_name = pass_name + "_pass"
    for name, fx in globals().items():
        if name == custom_pass_name:
            if callable(fx) and fx.__module__ == __name__:
                return fx
    return None

def compute_passes(args):
    passes = MAIN_PASSES.copy()
    flags = find_flags(args)
    for flag in AUX_PASSES.keys():
        if flag in flags: passes += [AUX_PASSES[flag]]
    # skip passes by name
    skipped_passes = get_flag_args('s', args)
    if skipped_passes:
        passes = handle_skip_flag(passes, skipped_passes)
    print("passes to run: ", passes)
    return passes

def execute_pass(pass_name, run_info):
    pass_function = get_custom_pass(pass_name)
    if pass_function: # if a custom pass function is defined, call it first
        pass_function(run_info)
    elif exists(f"scripts/{pass_name}.sh"): # otherwise look for a script
        call(["bash", f"scripts/{pass_name}.sh", run_info.run_name])
    else:
        print(f"I couldn't find a way to execute pass: {pass_name}\n"
                f"Either define a function called {pass_name}_pass in passes.py,\n"
                f"or create a bash script at scripts/{pass_name}.sh")

def execute_passes(run_path, args):
    run_info = RunInfo(run_path, args)
    for pass_name in compute_passes(args):
        execute_pass(pass_name, run_info)
