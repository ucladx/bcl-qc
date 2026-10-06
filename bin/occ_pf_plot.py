#!/usr/bin/env python3
# Plot InterOp occupancy against pass-filter percentage by lane.
import os, sys, logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_run_metrics(run_dir):
    from interop import py_interop_run_metrics, py_interop_run, py_interop_table
    from numpy import zeros, float32
    import pandas as pd
    run_metrics = py_interop_run_metrics.run_metrics()
    valid = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid[py_interop_run.ExtendedTile] = 1
    valid[py_interop_run.Tile] = 1
    valid[py_interop_run.Extraction] = 1
    try:
        run_metrics.read(run_dir, valid)
    except Exception as e:
        logging.error(f"Error occurred trying to open {run_dir}: {e}")
        return None
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    n = columns.size()
    if n == 0:
        logging.warning("No interop data found in run metrics.")
        return None
    headers = []
    for i in range(n):
        c = columns[i]
        if c.has_children():
            headers.extend([f"{c.name()} ({s})" for s in c.subcolumns()])
        else:
            headers.append(c.name())
    count = py_interop_table.count_table_columns(columns)
    offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, offsets)
    data = zeros((offsets.size(), count), dtype=float32)
    py_interop_table.populate_imaging_table_data(run_metrics, columns, offsets, data.ravel())
    return pd.DataFrame(data, columns=headers)

def main(run_dir, out_dir):
    df = parse_run_metrics(run_dir)
    if df is None:
        logging.warning("Unable to parse Interop files, skipping plot generation.")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.scatterplot(data=df, x="% Occupied", y="% Pass Filter", hue="Lane", alpha=0.5, s=8)
    plt.xlim([0, 100])
    plt.ylim([50, 100])
    plt.legend(title="Lane", bbox_to_anchor=[1.2, 0.9])
    plt.tight_layout()
    path = os.path.join(out_dir, "occ_pf_lane_mqc.jpg")
    logging.info(f"Saving Occupied vs Pass Filter graph to: {path}")
    plt.savefig(path, dpi=300)
    plt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
