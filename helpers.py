import os
import re
import sys
import pandas as pd
from subprocess import call
import matplotlib.pyplot as plt
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32

# This file includes helper functions used by bclqc.py

def get_run_path():
    """
    Returns the run path from the command line.
    Assumes usage: `python3 bclqc.py ... run_path`
    """
    run_path = sys.argv[-1]
    if run_path[-1] != '/':  run_path += '/'
    return run_path

def shell_exec(*args):
    call(["bash", *args])

def get_script(script_name):
    return f"scripts/{script_name}.sh"

def exec_pass(script_name, *args):
    shell_exec(get_script(script_name), *args)

def is_samplesheet(file_name):
    return re.search("^SampleSheet_.*csv$", file_name)

def get_index(file_name):
    # assumes samplesheet format is "SampleSheet_{index}.csv"
    return file_name[12:-4]

def get_indices(run_path):
    indices = [get_index(f) for f in os.listdir(run_path) if is_samplesheet(f)]
    if indices:
        return indices
    else:
        print("No samplesheets found in " + run_path)
        return None

def get_run_id(run_path: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2/',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    return run_path.split('/')[-2]

def get_sample_ids(fastq_list):
    fastq_list_df = pd.read_csv(fastq_list)
    return set(fastq_list_df['RGSM'])

# creates a dict mapping command line flags to their arguments
# e.x. {'O': ['a', 'b', 'c'], 'B': ['bed.BED']}
def parse_args(cli_args):
    args_dict = {}
    current_flag = ""
    for arg in cli_args:
        if arg.startswith('-'):
            current_flag = arg[1:]
            args_dict[current_flag] = []
        else:
            args_dict[current_flag] += [arg]
    return args_dict

def get_custom_passes():
    """
    Returns the list of custom passes from the command line.
    Assumes usage: `python3 bclqc.py -P pass1 pass2 ... run_path`
    """
    args = parse_args(sys.argv[1:-1])
    return args.get('P')

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

    return pd.DataFrame(data, columns=headers)

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

        save_dir = f"/mnt/pns/bams/{get_run_id(run_path)}/"
        image_path = save_dir + f"occ_pf_{view.lower()}_mqc.jpg"
        print("saving occ pf graph to " + image_path)
        plt.savefig(image_path, dpi=300)
        plt.close()
