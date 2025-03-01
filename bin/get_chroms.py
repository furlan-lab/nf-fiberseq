#!/usr/bin/env python
import re
import logging
import sys
import argparse
import os
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Compute coverage statistics")
    parser.add_argument("--input", required=True, help="Genome reference file")
    parser.add_argument("--min_contig_length", type=int,required=True, help="Get minimum contig length for filtering")
    return parser.parse_args()

args = get_args()

FIRST_REPORT = True
ref = args.input
def get_ref():   
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)


def get_fai():
    fai = f"{get_ref()}.fai"
    if not os.path.exists(fai):
        raise ValueError(f"FIRE: reference index file {fai} does not exist")
    return fai

def get_fai_df():
    fai = get_fai()
    return pd.read_csv(fai, sep="\t", names=["chr", "length", "x", "y", "z"])

FAI_DF = get_fai_df()

def get_chroms():
    global FIRST_REPORT
    min_contig_length = args.min_contig_length
    skipped_contigs = FAI_DF["chr"][FAI_DF["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    chroms = FAI_DF["chr"][FAI_DF["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    # chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {FAI_DF['chr']}"
        )
    return chroms

if __name__ == "__main__":
    chrom_list = get_chroms()
    print("\n".join(chrom_list))