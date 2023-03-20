import sys
from os.path import exists
from passes import execute_passes

"""
Usage: python3 bcl-qc.py run_path [flags...]

flags:
    -b ----- Specify a custom bed file for the alignment step
    -s ----- Skip passes by name or use "-s MAIN" to skip all main passes
    -M ----- Run MegaQC on this run
    -a ----- Upload the results of this run (BAMs, gVCFs) to Azure
"""


def samplesheet_missing(run_path, barcode_format):
    return not exists(run_path + f"SampleSheet_{barcode_format}.csv")

def generate_samplesheet(run_path, barcode_format):
    # TODO generate samplesheet based on corresponding Beaker extract
    # for now, error out appropriately
    raise Exception(f"SampleSheet_{barcode_format}.csv not found in {run_path}")

def qc_run(run_path: str, args: str):
    print(f"bcl-qc on {run_path}")
    # if samplesheet_missing(run_path, barcode_format): generate_samplesheet(run_path, barcode_format)
    execute_passes(run_path, args)

if __name__ == "__main__":
    # absolute path to a specific run's output directory
    # run_path = sys.argv[-1]
    # if run_path[-1] != '/':  run_path += '/'
    # args = sys.argv[1:-1]
    # qc_run(run_path, args)

### TODO ###
# Fix run watcher
# Better file organization
# Sample specific bed files (parse samplesheet)
# Support passes config file? 

    

# Unit tests:
#     sample sheet generation
#     run detection
#     Pass manipulation
