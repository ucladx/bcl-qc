import sys
from os.path import exists
from passes import exec_passes, get_barcode_format

def samplesheet_missing(run_path, barcode_format):
    return not exists(run_path + f"SampleSheet_{barcode_format}.csv")

def generate_samplesheet(run_path, barcode_format):
    # TODO generate samplesheet based on corresponding Beaker extract
    # for now, error out appropriately
    raise Exception(f"SampleSheet_{barcode_format}.csv not found in {run_path}")

def qc_run(run_path: str, flags: str):
    barcode_format = get_barcode_format(flags)
    print(f"bcl-qc on {run_path}:\nflags: {flags}\nbarcode: {barcode_format}")
    if samplesheet_missing(run_path, barcode_format): generate_samplesheet(run_path, barcode_format)
    exec_passes(run_path, flags)

if __name__ == "__main__":
    # absolute path to a specific run's output directory
    run_path = sys.argv[-1]
    if run_path[-1] != '/':  run_path += '/'

    flags = sys.argv[1:-1]
    qc_run(run_path, flags)
