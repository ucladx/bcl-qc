import sys
from os.path import exists
from passes import exec_pass, compute_passes, RunInfo

def samplesheet_missing(run_info):
    run_path = run_info.run_path
    if run_path[-1] != '/': run_path += '/' # ensure trailing slash 
    return not exists(run_path + f"SampleSheet_{run_info.barcode}.csv")

def generate_samplesheet(run_info):
    # TODO generate samplesheet based on corresponding Beaker extract
    # for now, error out appropriately
    raise Exception(f"SampleSheet_{run_info.barcode}.csv not found in {run_path}")

def qc_run(run_path: str, flags: str):
    run_info = RunInfo(flags, run_path)
    print(f"bcl-qc on {run_path}:\nflags: {flags}\nbarcode: {run_info.barcode}")
    if samplesheet_missing(run_info): generate_samplesheet(run_info)
    for pass_name in compute_passes(flags):
        exec_pass(pass_name, run_info)

if __name__ == "__main__":
    # ex: python3 bcl-qc.py -u -m ~/221013_A01718_0014_AHNYGGDRX2
    run_path = sys.argv[-1]
    flags = sys.argv[1:-1]
    qc_run(run_path, flags)
