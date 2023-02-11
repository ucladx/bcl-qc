# Main passes run by default
MAIN_PASSES = [
    "demux", # bash
    "align", # bash
    ]

AUX_PASSES = {
    'm' : "multiqc", # bash
    'M' : "megaqc", # bash
    'a' : "azure", # py
}

class RunInfo:
    def __init__(self, flags, run_path):
        self.barcode = get_barcode_format(flags)
        self.run_path = run_path
        self.run_name = get_run_name(run_path)

def azure_pass(run_info):
    print("Azure Upload\n-------------------------")
    # azure upload logic goes here
    return

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
    save_occ_pf_plot(run_info.run_path)
    call(["bash", "~/scripts/multiqc.sh", run_info.barcode, run_info.run_name])

def get_flag_args(char, flags):
    flag = "-" + char
    if flag in flags:
        start = flags.index(flag) + 1
        for ele in flags[start:]:
            if ele.startswith('-'):
                end = flags.index(ele) - 1
                return [arg for arg in flags[start:end]]
    return []

def compute_passes(flags):
    passes = MAIN_PASSES
    all_flags = ''.join(flags[0:-1])
    for aux_pass_flag in AUX_PASSES.keys():
        if aux_pass_flag in all_flags: passes += AUX_PASSES[aux_pass_flag]
    if 's' in all_flags:
        for arg in get_flag_args('s', flags):
            if arg == "all":
                passes = []
            else if arg == "main":
                for pass_name in MAIN_PASSES:
                    passes.remove(pass_name)
            else:
                passes.remove(arg)
    return passes

def call_pass(pass_name, run_info):
    pass_function = globals()[pass_name + "_pass"] # scary --- we pull this function out of the global namespace
    pass_function(run_info.barcode, run_info.run_path)
