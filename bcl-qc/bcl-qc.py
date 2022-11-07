import matplotlib.pyplot as plt
import sys
from subprocess import Popen, PIPE
from numpy import zeros, float32
from pandas import DataFrame
from seaborn import scatterplot
from interop import py_interop_run_metrics, py_interop_run, py_interop_table, py_interop_summary

def parse_run_metrics(run_dir: str):
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
    run_metrics.read(run_dir, valid_to_load)
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
        plt.ylim([0, 100])
        plt.legend(title=hue, bbox_to_anchor=[1.2, 0.9])
        plt.tight_layout()
        plt.savefig(run_dir, "/occupancy" + f"_{hue.lower()}_mcq.jpg", dpi=600)
        plt.close()


def occ_pf(run_dir: str):
    """
    Saves a scatter plot of % Occupied x % Pass Filter by lane for the run in `run_dir`
    """
    df = parse_run_metrics(run_dir)
    if df:
        save_occ_pf_plot(df)
    else:
        raise Exception("No data available from given InterOp files.")

# TODO comments from readme
def demux(run_dir: str):
    Popen(["perl",
           "-a",
           "-F'\\t'",
           "-ne",
           "'BEGIN{%m=map{@c=split(\"\\t\"); (\"$c[0] $c[2]\", \"$c[3],$c[4]\")}`cat ~/illumina_barcodes.txt`} chomp($F[1]); print join(\",\",$F[0],$m{$F[1]})'",
           "sample_index_map.txt"])
    Popen(["rsync",
           "-rhlW",
           "--info=progress2",
           "--exclude 'Thumbnail_Images' /mnt/pns/runs/{} /staging/hot/"])
    Popen(["mkdir",
           f"/staging/hot/reads/{run_dir}"])
    Popen(["dragen",
           "--bcl-conversion-only true",
           "--bcl-use-hw false",
           "--bcl-only-matched-reads true",
           f"--bcl-input-directory /staging/hot/{run_dir}",
           f"--sample-sheet /staging/hot/{run_dir}/SampleSheet_I10.csv ",
           f"--output-directory /staging/hot/reads/{run_dir}/I10"])

# TODO comments from readme
def align(run_dir: str):
    run_dir = run_dir[1:] # strip leading '/'
    BED = "/mnt/pns/tracks/ucla_mdl_cancer_ngs_v1_exon_targets.hg38.bed"
    fastq_list = f"/staging/hot/reads/{run_dir}/I10/Reports/fastq_list.csv"
    bams_dir = f"/mnt/pns/bams/{run_dir}/{{}}"

    cut1 = Popen(["cut", "-f2", "-d", fastq_list], stdout=PIPE)
    grep1 = Popen(["grep", "-Ev", "^RGSM"], stdin=cut1.stdout, stdout=PIPE)
    sort1 = Popen(["sort", "-u"], stdin=grep1.stdout, stdout=PIPE)
    xargs1 = Popen(["xargs", "-L1", "-I{}", f"mkdir -p {bams_dir}"], stdin=sort1.stdout)

    cut2 = Popen(["cut", "-f2", "-d", fastq_list], stdout=PIPE)
    grep2 = Popen(["grep", "-Ev", "^RGSM"], stdin=cut2.stdout, stdout=PIPE)
    sort2 = Popen(["sort", "-u"], stdin=grep2.stdout, stdout=PIPE)
    xargs2 = Popen(["xargs", "-L1", "-I{}",
                    ("dragen --intermediate-results-dir /staging/tmp"
                            "--enable-map-align true"
                            "--enable-map-align-output true"
                            "--output-format BAM"
                            "--enable-duplicate-marking true"
                            "--generate-sa-tags true"
                            "--enable-sort true"
                            "--ref-dir /staging/human/reference/hg38_alt_masked_graph_v2"
                            "--qc-coverage-tag-1 exon"
                            f"--qc-coverage-region-1 {BED}"
                            "--qc-coverage-reports-1 cov_report"
                            "--qc-coverage-ignore-overlaps true"
                            "--enable-variant-caller true"
                            "--vc-emit-ref-confidence GVCF"
                            "--enable-hla true"
                            f"--fastq-list {fastq_list}"
                            "--fastq-list-sample-id {}"
                            f"--output-directory {bams_dir} --output-file-prefix {{}}")]
                    ,stdin=sort2.stdout)

def multiqc_report(run_dir: str):
    run_dir = run_dir[1:] # strip leading '/'
    Popen(["rm", "-f", f"/mnt/pns/bams/{run_dir}/*/*.wgs_*.csv"])
    Popen(["multiqc",
           f"--outdir /mnt/pns/bams/{rundir}",
           f"/mnt/pns/bams/{run_dir}",
           f"/staging/hot/reads/{run_dir}/I10"])


def qc_run(run_dir: str):
    passes = [occ_pf, demux, align, multiqc_report]
    for qc_pass in passes:
        qc_pass(run_dir)

if __name__ == "__main__":
    run_dir = sys.argv[1]
    qc_run(run_dir)

