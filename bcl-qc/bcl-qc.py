import matplotlib.pyplot as plt
import sys
from subprocess import call, PIPE
from numpy import zeros, float32
from pandas import DataFrame
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table

def parse_run_metrics(run_path: str):
    """
    Returns a dataframe of `interop_imaging_table` or returns None if no data is found

    More info on the `interop` library:
    http://illumina.github.io/interop/index.html
    """
    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    # summary = py_interop_summary.run_summary()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    run_metrics.read(run_path, valid_to_load)
    # py_interop_summary.summarize_run_metrics(run_metrics, summary)

    # Create the columns
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)

    headers = []
    for i in range(columns.size()):
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
    py_interop_table.populate_imaging_table_data(
        run_metrics, columns, row_offsets, data.ravel()
    )

    # Make a DataFrame
    df = DataFrame(data, columns=headers)

    # Return None if there is no data
    if df.shape[0] == 0:
        return None
    else:
        return df

# Given an interop imaging table dataframe,
# saves a scatter plot of its Occupied% x PassFilter% data to `run_dir`
def save_occ_pf_plot(df, run_dir: str):
    x = "% Pass Filter"
    y = "% Occupied"
    hues = ["Lane"] # can add e.g. "Tile" or "Cycle" for different views

    for hue in hues:
        scatterplot(
            data=df,
            x=x,
            y=y,
            hue=hue,
            alpha=0.5,
            s=8,
        )
        plt.xlim([0, 100])
        plt.ylim([50, 100])
        plt.legend(title=hue, bbox_to_anchor=[1.2, 0.9])
        plt.tight_layout()
        print("saving occ pf graph")
        plt.savefig(f"/staging/hot/reads/{run_dir}" + "occupancy" + f"_{hue.lower()}_mqc.jpg", dpi=600)
        plt.close()

def occ_pf(run_path: str):
    """
    Saves a scatter plot of % Occupied x % Pass Filter by lane for the run in `run_dir`
    """
    df = parse_run_metrics(run_path)
    if df is not None:
        run_dir = dir_from_path(run_path)
        save_occ_pf_plot(df, run_dir)
    else:
        raise Exception("No data available from given InterOp files.")

def dir_from_path(path: str):
    # strip everything but dir name
    return [x for x in path.split('/') if x][-1]

def qc_run(run_path: str):
    occ_pf(run_path)
    run_dir = dir_from_path(run_path)
    call(["bash", "bcl-qc.sh", run_dir])

if __name__ == "__main__":
    run_path = sys.argv[1]
    qc_run(run_path)

