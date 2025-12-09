// modules/interop_plot.nf - Generate Occupied vs Pass Filter plot from Illumina InterOp

process INTEROP_PLOT {
    tag "${run_dir.name}"
    label 'low_memory'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path run_dir

    output:
    path "occ_pf_lane_mqc.jpg", emit: plot, optional: true

    script:
    """
    #!/usr/bin/env python3
    import os
    import sys

    try:
        from interop import py_interop_run_metrics, py_interop_run, py_interop_table
        import pandas as pd
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as e:
        print(f"Warning: Could not import required module: {e}", file=sys.stderr)
        sys.exit(0)  # Exit gracefully - plot is optional

    def parse_run_metrics(run_dir):
        run_metrics = py_interop_run_metrics.run_metrics()
        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        valid_to_load[py_interop_run.ExtendedTile] = 1
        valid_to_load[py_interop_run.Tile] = 1
        valid_to_load[py_interop_run.Extraction] = 1

        try:
            run_metrics.read(run_dir, valid_to_load)
        except Exception as e:
            print(f"Warning: Could not read InterOp data: {e}", file=sys.stderr)
            return None

        columns = py_interop_table.imaging_column_vector()
        py_interop_table.create_imaging_table_columns(run_metrics, columns)
        ncolumns = columns.size()

        if ncolumns == 0:
            return None

        headers = []
        for i in range(ncolumns):
            column = columns[i]
            if column.has_children():
                headers.extend([f"{column.name()} ({subname})" for subname in column.subcolumns()])
            else:
                headers.append(column.name())

        column_count = py_interop_table.count_table_columns(columns)
        row_offsets = py_interop_table.map_id_offset()
        py_interop_table.count_table_rows(run_metrics, row_offsets)
        data = np.zeros((row_offsets.size(), column_count), dtype=np.float32)

        py_interop_table.populate_imaging_table_data(run_metrics, columns, row_offsets, data.ravel())

        return pd.DataFrame(data, columns=headers)

    # Main
    df = parse_run_metrics('${run_dir}')

    if df is None or df.empty:
        print("No InterOp data available, skipping plot", file=sys.stderr)
        sys.exit(0)

    x = "% Occupied"
    y = "% Pass Filter"

    if x not in df.columns or y not in df.columns:
        print(f"Required columns not found in InterOp data", file=sys.stderr)
        sys.exit(0)

    sns.scatterplot(data=df, x=x, y=y, hue="Lane", alpha=0.5, s=8)
    plt.xlim([0, 100])
    plt.ylim([50, 100])
    plt.legend(title="Lane", bbox_to_anchor=[1.2, 0.9])
    plt.tight_layout()
    plt.savefig("occ_pf_lane_mqc.jpg", dpi=300)
    plt.close()

    print("Generated occ_pf_lane_mqc.jpg")
    """
}
