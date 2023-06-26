import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import call
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32
from bclqc import MAIN_PASSES

# This file includes helper functions used by bclqc.py

def trim_slash(path):
    if path[-1] == '/': path = path[:-1]
    return path

def get_run_name(run_dir: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    split_fields = run_dir.split('/')
    if split_fields[-1] == '':
        return split_fields[-2]
    else:
        return split_fields[-1]

def get_arg_dirs():
    if len(sys.argv) < 4:
        raise Exception("Usage: python3 bclqc.py [options] fastqs_dir bams_dir run_dir")
    run_dir = sys.argv[-1]
    bams_dir = sys.argv[-2]
    fastqs_dir = sys.argv[-3]
    run_name = get_run_name(run_dir)
    dirs = [trim_slash(d) for d in [fastqs_dir, bams_dir, run_dir]]
    for i, dir in enumerate(dirs):
        if run_name not in dir:
            dirs[i] = f"{dir}/{run_name}"
    return dirs

def get_cli_args():
    return parse_args(sys.argv[1:-3])

def get_custom_passes():
    return get_cli_args().get('P')

def get_passes():
    custom_passes = get_custom_passes()
    return custom_passes if custom_passes else MAIN_PASSES

def print_cmd(cmd):
    print(" ".join(cmd) + "\n")

def get_exec_cmd():
    if '-d' in sys.argv or '--dry-run' in sys.argv:
        return print_cmd
    else:
        return call

def exec_pass(script_name, *args):
    call(["bash", f"scripts/{script_name}.sh", *args])

def is_samplesheet(file_name):
    return re.search("^SampleSheet_.*csv$", file_name)

def get_index(file_name):
    # assumes samplesheet format is "SampleSheet_{index}.csv"
    return file_name[12:-4]

def get_indices(run_dir):
    indices = [get_index(f) for f in os.listdir(run_dir) if is_samplesheet(f)]
    if not indices:
        raise("No samplesheets found in " + run_dir)
    return indices

def get_sample_ids(fastq_list):
    fastq_list_df = pd.read_csv(fastq_list)
    return set(fastq_list_df['RGSM'])

# creates a dict mapping command line flags to their arguments
# e.x. {'P': ['demux', 'align']}
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
    Assumes usage: `python3 bclqc.py -P pass1 pass2 ... fastqs_dir bams_dir run_dir`
    """
    args = parse_args(sys.argv[1:-3])
    return args.get('P')

def parse_run_metrics(run_dir: str):
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
        run_metrics.read(run_dir, valid_to_load)
    except:
        print(f"Error occured trying to open {run_dir}")
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

def save_occ_pf_plot(run_dir, output_dir, exec_cmd=call):
    """
    Saves a % Occupied x % Pass Filter scatter to `SAVE_DIR`.
    """
    df = parse_run_metrics(run_dir)
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

        image_path = output_dir + f"/occ_pf_{view.lower()}_mqc.jpg"
        if exec_cmd == call:
            print("saving occ pf graph to " + image_path)
            plt.savefig(image_path, dpi=300)
        else:
            print("DRY RUN: saving occ pf graph to " + image_path)
        plt.close()
