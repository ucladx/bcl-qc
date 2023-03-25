#!/usr/bin/env python3

import pandas as pd
import csv
from subprocess import run
from Bio.Seq import Seq
from skbio.sequence import DNA

# TODO: Don't use Perl for this
# Instead, call `samplesheet.py` from `new_run_watcher.sh` if SampleSheet.xml is not in the run folder
# If we find that no SampleSheet already exists, then run the equivalent of this command (in Python)

# def sample_sheet_cmd():
#     Popen(["perl",
#         "-a",
#         "-F'\\t'",
#         "-ne",
#         "'BEGIN{%m=map{@c=split(\"\\t\"); (\"$c[0] $c[2]\", \"$c[3],$c[4]\")}`cat ~/illumina_barcodes.txt`} chomp($F[1]); print join(\",\",$F[0],$m{$F[1]})'",
#         "sample_index_map.txt"])

def get_dragen_version():
    """
    Attempts to get the dragen version on the current machine
    Returns None if a problem occurs
    """
    try:
        # run dragen and get the relevant output from `CompletedProcess`
        version_num = run(["dragen", "--version"], capture_output=True).stdout.split()[-1]
    except:
        print("Problem running \'dragen --version\' or parsing its output")
        return None
    
    # given b'07.021.645.4.0.3', return '4.0.3'
    version_num = version_num.split(b".", 3)[-1].decode("utf-8")
    return version_num

def create_samplesheet_header():
    """
    Creates the header for the BCL-QC SampleSheet
    """
    # TODO handle SampleSheets that aren't I10 (add OverrideCycles to bclconvert_string)
    header_string = "[Header]\nFileFormatVersion,2\n"
    reads_string = "[Reads]\nRead1Cycles,151\nRead2Cycles,151\nIndex1Cycles,10\nIndex2Cycles\n"
    bclconvert_string = f"[BCLConvert_Settings]\nSoftwareVersion,{get_dragen_version()}\nBarcodeMismatchesIndex1,1\nBarcodeMismatchesIndex2,1\n"
    return header_string + reads_string + bclconvert_string

def beaker_to_samplesheet(beaker_df: pd.DataFrame):
    """
    Given a DataFrame of the Beaker extract,
    returns a string representing the BCL-QC SampleSheet CSV
    """

    samplesheet_cols = [
        "sample_id",
        "index1",
        "index2",]

    samplesheet_df = pd.DataFrame(columns=samplesheet_cols)
    samplesheet_df['sample_id'] = beaker_df['sample_id']
    samplesheet_df['index1', 'index2'] = get_indices_pairs(beaker_df[['plate_id', 'sample_id']])
    samplesheet = create_samplesheet_header() + samplesheet_df
    return samplesheet

def beaker_to_minimal_input(beaker_df: pd.DataFrame):
    min_input_cols = [
        "sample_id",
        "plate_id",
        "well_id",
        "run_id",
        "specimen_type",
        "tissue_source",
        "patient_name",
        "sex",
        "MRN",
        "dob",
        "family_id",
        "provider_name",
        "test_ordered",
        "diagnosis",
        "pre_hyb conc",
        "pre_hyb_avg_size",
        "final_lib_conc",
        "final_lib_size",
    ]

    min_input_df = pd.DataFrame(columns=min_input_cols)
    # copy columns that have the same names in both DFs
    for dup_col_name in ['plate_id', 'well_id']:
        min_input_df[f'{dup_col_name}'] = beaker_df[f'{dup_col_name}']

    # TODO manually map the columns that don't have the same name
    
    raise NotImplementedError

def get_indices_pairs(beaker_df: pd.DataFrame):
    """
    Given the Beaker DataFrame, returns a list of (index, index2),
    corresponding to each row
    """
    BARCODES_PATH = '../static/illumina_barcodes.txt'
    try:
        ib_df = pd.read_csv(BARCODES_PATH,
                            sep='\t',
                            header=None,
                            colums=['plate_id', 'index_id', 'sample_id', 'index1', 'index2'])
        indices = []
        for beaker_row in beaker_df.itertuples():
            for ib_row in ib_df.itertuples():
                if beaker_row.plate_id == ib_row.plate_id and beaker_row.sample_id == ib_row.sample_id:
                    indices.append((ib_row.index1, ib_row.index2))
    except:
        print(f"Failed to open {BARCODES_PATH}")

def create_samplesheet(beaker_path: str):
    try:
        beaker_df = pd.read_csv(beaker_path, sep = '\t')
        if beaker_df.empty:
            print(f"No data found from Beaker extract at {beaker_path}")
            return
        else:
            return beaker_to_samplesheet(beaker_df)
    except:
        print(f"Error opening {beaker_path}")

def save_samplesheet(run_path: str, beaker_path: str):
    """
    Generates the BCL-QC SampleSheet
    and saves it to `run_path/SampleSheet_I10.csv`
    """
    sample_sheet = create_samplesheet(beaker_path)
    try:
        file = open(f"{run_path}/SampleSheet_I10.csv", "w")
        file.write(sample_sheet)
        file.close()
    except:
        print(f"Failed to write to {run_path}/SampleSheet_I10.csv")
        return

def check_barcodes():
    """
    Check for possible collisions between given barcodes
    """
    pass

# utilities for identifying barcode similarity/collision
def reverse_complement(index):
    rc = Seq(index).reverse_complement()
    return(rc)

def to_rc(index_series):
    rc = index_series.apply(lambda x : Seq(x).reverse_complement(),axis=1)
    return(rc)

def get_hamming(a,b):
    a = DNA(a)
    hamming_distance = a.distance(b)
    return(hamming_distance)

def get_mismatches(a,b):
    hamming_distance = get_hamming(a,b)
    mismatches = round(len(a)*hamming_distance)
    return(mismatches)

def mismatch_check(a,b,n=3):
    mismatches = get_mismatches(a,b)
    passed = n > mismatches
    return(passed)

def mismatch_vector(a,series_b):
    mismatches = series_b.apply(lambda x : get_mismatches(a,x))
    mismatches = pd.Series(mismatches)
    mismatches.index = series_b
    mismatches.name = a
    return(mismatches)

def pairwise_mismatch_check(series_a,series_b):
    mismatches = series_a.apply(lambda x : mismatch_vector(x,series_b))
    mismatches.index = series_a.values
    #mismatches.index = mismatches.columns
    return(mismatches)

#debug utils
def fetch_ib():
    out = pd.read_csv("illumina_barcodes.txt",sep="\t")
    return(out)

def combined_barcode(df):
    out = df["index"] + df["index2"]
    return(out)

def add_combined_barcode(df):
    df["combined"] = combined_barcode(df) 
    return(df)

def ib():
    out = fetch_ib()
    out = add_combined_barcode(out)
    return(out)

def archer():
    out = ib()
    out = out[out.Plate_ID.str.contains("FusionPlex")]
    return(out)

def twist():
    out = ib()
    out = out[out.Plate_ID.str.contains("Twist")]
    return(out)

def main():
    a = archer()
    t = twist()
    amm = pairwise_mismatch_check(a.combined,a.combined)
    tmm = pairwise_mismatch_check(t.combined,t.combined)
    amm.to_csv("archer_mms.txt")
    tmm.to_csv("twist_mms.txt")

#dump mismatch stats 
if __name__ == "__main__":
    main()
