import matplotlib.pyplot as plt
import sys
from subprocess import call, PIPE
from numpy import zeros, float32
from pandas import DataFrame
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from os.path import exists

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

def occ_pf_plot(df: DataFrame, run_path: str):
    """
    Given an interop imaging table dataframe,
    saves a % Occupied x % Pass Filter scatter to `SAVE_DIR`.
    """
    run_name = get_run_name(run_path)
    SAVE_DIR = f"/staging/hot/reads/{run_name}/I10/"

    x = "% Pass Filter"
    y = "% Occupied"
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

        image_path = SAVE_DIR + f"occ_pf_{view.lower()}_mqc.jpg"
        print("saving occ pf graph to " + image_path)
        plt.savefig(image_path, dpi=300)
        plt.close()

def get_run_name(run_path: str):
    """
    Given '.../221013_A01718_0014_AHNYGGDRX2/',
    returns '221013_A01718_0014_AHNYGGDRX2'
    """
    return [x for x in run_path.split('/') if x][-1]

def samplesheet_exists(run_path):
    if run_path[-1] != '/': run_path += '/' # ensure trailing slash 
    return exists(run_path + "SampleSheet_I10.csv")

def qc_run(run_path: str):
    # check for 
    if samplesheet_exists(run_path):
        df = parse_run_metrics(run_path)
        if df is not None:
            occ_pf_plot(df, run_path)
        else:
            print("Unable to parse Interop files ---",
                "could not generate % Occupied x % Pass Filter graph.")

        # ex: `bash bcl-qc.sh 221013_A01718_0014_AHNYGGDRX2`
        call(["bash", "bcl-qc.sh", get_run_name(run_path)])
    else:
        print(f"SampleSheet_I10.csv not found in {run_path}")

if __name__ == "__main__":
    run_path = sys.argv[1]
    qc_run(run_path)
