import matplotlib.pyplot as plt
import os
import re
from pandas import DataFrame
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32
from subprocess import call
from os.path import exists

# Main passes are run by default
MAIN_PASSES = [
    "demux",
    "align",
]

# Aux passes are defined with a flag that can be given to enable them
AUX_PASSES = {
    'm' : "multiqc",
    'M' : "megaqc",
    'a' : "azure",
}
AUX_PASS_FLAGS = AUX_PASSES.keys()

# TODO just get all barcode formats based on SampleSheet_IDX existing in run folder
DEFAULT_BARCODE_FORMAT = "I10" # default sample sheet barcode_format, e.g. I10, I8N2
DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"


def is_samplesheet(file_name):
    return re.search("^SampleSheet.*", file_name)

def get_barcode(file_name):
    # assumes samplesheet format is "SampleSheet_{index}.csv"
    return file_name[12:-4]

def get_barcodes(run_path):
    return [get_barcode(f) for f in os.listdir(run_path) if is_samplesheet(f)]

def get_barcode_format(args):
    barcode_format = DEFAULT_BARCODE_FORMAT
    custom_format = get_flag_args('i', args)
    if custom_format: barcode_format = custom_format[0]
    return barcode_format

def get_run_name(run_path: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2/',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    return [x for x in run_path.split('/') if x][-1]

class RunInfo:
    def __init__(self, run_path, args):
        self.run_path = run_path
        self.args = args
        self.barcodes = get_barcode_format(args)
        self.run_name = get_run_name(run_path)

class PassInfo:
    def __init__(self, args, pass_name, pass_function):
        self.args = args
        self.pass_name = pass_name
        self.pass_function = pass_function

def azure_pass(run_info):
    print("TEST Azure Upload\n-------------------------")
    return

def align_pass(run_info):
    bed_path = DEFAULT_BED_PATH
    user_bed_path = get_flag_args('b', run_info.args)
    if user_bed_path: bed_path = user_bed_path[0]
    call(["bash", "/home/iatol/bclqc/bcl-qc/scripts/align.sh", run_info.barcode, run_info.run_name, bed_path])

def parse_run_metrics(run_path: str):
    """
    Returns a DataFrame of `interop_imaging_table` or returns None if no data is found.

    More info on `interop_imaging_table` and the `interop` library:
    http://illumina.github.io/interop/index.html
    https://github.com/Illumina/interop
    """

    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    try:
        run_metrics.read(run_path, valid_to_load)
    except:
        print(f"Error occured trying to open {run_path}")
        return None

    # Set up data table
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    ncolumns = columns.size()
    if ncolumns == 0: return None # no data

    headers = []
    for i in range(ncolumns):
        column = columns[i]
        if column.has_children():
            headers.extend(
                [f"{column.name()} ({subname})" for subname in column.subcolumns()])
        else:
            headers.append(column.name())

    column_count = py_interop_table.count_table_columns(columns)
    row_offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, row_offsets)
    data = zeros((row_offsets.size(), column_count), dtype=float32)

    # Populate table
    py_interop_table.populate_imaging_table_data(
        run_metrics, columns, row_offsets, data.ravel()
    )

    return DataFrame(data, columns=headers)

def save_occ_pf_plot(run_path: str):
    """
    Saves a % Occupied x % Pass Filter scatter to `SAVE_DIR`.
    """

    df = parse_run_metrics(run_path)
    if df is None:
        print("Unable to parse Interop files.")
        return

    x = "% Occupied"
    y = "% Pass Filter"
    views = ["Lane"] # can add "Tile" or "Cycle" for more views
    for view in views:
        scatterplot(
            data=df,
            x=x,
            y=y,
            hue=view,
            alpha=0.5,
            s=8,
        )
        plt.xlim([0, 100])
        plt.ylim([50, 100])
        plt.legend(title=view, bbox_to_anchor=[1.2, 0.9])
        plt.tight_layout()

        save_dir = f"/mnt/pns/bams/{get_run_name(run_path)}/"
        image_path = save_dir + f"occ_pf_{view.lower()}_mqc.jpg"
        print("saving occ pf graph to " + image_path)
        plt.savefig(image_path, dpi=300)
        plt.close()

def multiqc_pass(run_info):
    print("MultiQC\n-------------------------")
    save_occ_pf_plot(run_info.run_path)
    call(["bash", f"/home/iatol/bclqc/bcl-qc/scripts/multiqc.sh", run_info.run_name])

def is_flag(arg):
    return arg.startswith('-')

def get_flag_args(char, args):
    flag = "-" + char
    if flag in args:
        start = args.index(flag) + 1
        if len(args[start:]) == 1:
            return args[start:]
        else:
            for arg in args[start:]:
                if is_flag(arg):
                    end = args.index(arg)
                    return args[start:end]
                elif arg == args[-1]:
                    end = args.index(arg) + 1
                    return args[start:end]
    return []

def find_flags(args):
    return [arg[1] for arg in args if is_flag(arg)]

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

def compute_passes(args):
    passes = MAIN_PASSES.copy()
    flags = find_flags(args)
    for flag in AUX_PASS_FLAGS:
        if flag in flags: passes += [AUX_PASSES[flag]]
    # skip passes by name
    skipped_passes = get_flag_args('s', args)
    if skipped_passes:
        passes = handle_skip_flag(passes, skipped_passes)
    print("passes to run: ", passes)
    return passes

def get_custom_pass(pass_name):
    custom_pass_name = pass_name + "_pass"
    for name, fx in globals().items():
        if name == custom_pass_name:
            if callable(fx) and fx.__module__ == __name__:
                return fx
    return None

def execute(pass_name, run_info):
    pass_function = get_custom_pass(pass_name)
    if pass_function: # if a custom pass function is defined, call it first
        pass_function(run_info)
    elif exists(f"/home/iatol/bclqc/bcl-qc/scripts/{pass_name}.sh"): # otherwise look for a script
        call(["bash", f"/home/iatol/bclqc/bcl-qc/scripts/{pass_name}.sh", run_info.barcode, run_info.run_name])
    else:
        print(f"I couldn't find a way to execute pass: {pass_name}\n"
                f"Either define a function called {pass_name}_pass in passes.py,\n"
                f"or create a bash script at /home/iatol/bclqc/bcl-qc/scripts/{pass_name}.sh")

def execute_passes(run_path, args):
    run_info = RunInfo(run_path, args)
    for pass_name in compute_passes(args):
        print(f"Running BCL-QC pass: {pass_name}")
        execute(pass_name, run_info)
