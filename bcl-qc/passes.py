import matplotlib.pyplot as plt
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

DEFAULT_BARCODE_FORMAT = "I10" # default sample sheet barcode_format, e.g. I10, I8N2
DEFAULT_BED_PATH = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"

def get_barcode_format(flags):
    barcode_format = DEFAULT_BARCODE_FORMAT
    custom_format = get_flag_args('i', flags)
    if custom_format: barcode_format = custom_format[0]
    return barcode_format

def get_run_name(run_path: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2/',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    return [x for x in run_path.split('/') if x][-1]

class RunInfo:
    def __init__(self, run_path, flags):
        self.run_path = run_path
        self.flags = flags
        self.barcode = get_barcode_format(flags)
        self.run_name = get_run_name(run_path)

def azure_pass(run_info):
    print("TEST Azure Upload\n-------------------------")
    return

def align_pass(run_info):
    bed_path = DEFAULT_BED_PATH
    user_bed_path = get_flag_args('b', run_info.flags)
    if user_bed_path: bed_path = user_bed_path
    call(["bash", "~/scripts/align.sh", run_info.barcode, run_info.run_name, bed_path])

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

        save_dir = f"/mnt/pns/runs/{get_run_name(run_path)}/"
        image_path = save_dir + f"occ_pf_{view.lower()}_mqc.jpg"
        print("saving occ pf graph to " + image_path)
        plt.savefig(image_path, dpi=300)
        plt.close()

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
                end = flags.index(ele)
                return [arg for arg in flags[start:end]]
    return []

def handle_skip_flag(passes, flags):
    new_passes = passes
    for arg in get_flag_args('s', flags):
        if arg == "all":
            return []
        elif arg == "main":
            for pass_name in MAIN_PASSES: new_passes.remove(pass_name)
        else:
            new_passes.remove(arg)
    return new_passes

def compute_passes(flags):
    passes = MAIN_PASSES
    found_flags = ''.join(flags)
    for flag in AUX_PASS_FLAGS:
        if flag in found_flags: passes += AUX_PASSES[flag]
    if 's' in found_flags: # skip passes by name
        passes = handle_skip_flag(passes, flags)
    return passes

def exec_pass(pass_name, run_info):
    if pass_name + "_pass" in dir(): # if a custom pass function is defined, call it first
        pass_function = globals()[pass_name + "_pass"] # scary --- we pull this function out of the global namespace
        pass_function(run_info)
    elif exists("~/scripts/{pass_name}.sh"): # otherwise look for a script
        call(["bash", "~/scripts/{pass_name}.sh", run_info.barcode, run_info.run_name])
    else:
        print(f"I couldn't find a way to execute pass: {pass_name}.\n"
                "Either define a function called {pass_name}_pass in passes.py,\n"
                "or create a bash script called {pass_name}.sh in ~/scripts/")

def exec_passes(run_path, flags):
    run_info = RunInfo(flags, run_path)
    for pass_name in compute_passes(flags):
        exec_pass(pass_name, run_info)