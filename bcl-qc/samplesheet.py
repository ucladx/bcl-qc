#!/usr/bin/env python3
import pandas as pd 
from Bio.Seq import Seq
from skbio.sequence import DNA

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
    mismatches.index = mismatches.columns
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
    out = add_combined_barcode(df)
    return(out)

def archer():
    out = ib()
    out = out[out.Plate_ID.str.contains("FusionPlex")]
    return(out)

def twist():
    out = ib()
    out = out[out.Plate_ID.str.contains("Twist")]
    return(out)
