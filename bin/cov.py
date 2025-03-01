#!/usr/bin/env python
import pandas as pd
import argparse
import sys
import math
import polars as pl
import re
import logging


def get_args():
    parser = argparse.ArgumentParser(description="Compute coverage statistics")
    parser.add_argument("--input", required=True, help="Path to bedgraph input file")
    parser.add_argument("--output_median", required=True, help="Output file directory for median coverage")
    parser.add_argument("--output_min", required=True, help="Output file directory for minimum coverage")
    parser.add_argument("--output_max", required=True, help="Output file directory for maximum coverage")
    parser.add_argument("--coverage_within_n_sd", type=float, required=True, help="Standard deviation multiplier")
    parser.add_argument("--min_coverage", type=float, required=True, help="Minimum allowed coverage")
    parser.add_argument("--chromosomes", required=True, help="Txt file with chromosomes list")
    return parser.parse_args()

args = get_args()

def read_chromosomes():
    with open(args.chromosomes, "r") as f:
        chromosomes = [line.strip() for line in f if line.strip()]  # Remove empty lines
    return chromosomes

chroms = read_chromosomes()

def get_min_coverage(median):
    sd = math.sqrt(median)
    mmin = median - args.coverage_within_n_sd * sd
    return max(mmin, args.min_coverage)


def get_max_coverage(median):
    sd = math.sqrt(median)
    return median + args.coverage_within_n_sd * sd


def weighted_median(df, val, weight):
    # group by value and sum the weights
    gdf = df.groupby(val)[weight].sum().reset_index().sort_values(val)
    print(gdf, file=sys.stderr)

    gdf["cumsum"] = gdf[weight].cumsum()
    gdf["cutoff"] = gdf[weight].sum() / 2.0
    print(gdf, file=sys.stderr)
    comparison = gdf[gdf["cumsum"] >= gdf["cutoff"]][val]
    # print(comparison, file=sys.stderr)
    return comparison.iloc[0]


def pandas_read():
    df = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "coverage"],
    )
    df = df.loc[df["coverage"] > 0]
    df = df.loc[df["chr"].isin(chroms)]
    df["weight"] = df["end"] - df["start"]
    return df


def polars_read():
    # Reading the CSV file using the lazy API
    df = (
        pl.read_csv(
            args.input,
            separator="\t",
            has_header=False,
            new_columns=["chr", "start", "end", "coverage"],
            low_memory=True,
        )
        .lazy()
        .filter(pl.col("coverage") > 0)
        .filter(pl.col("chr").is_in(chroms))
        .drop("chr")
        .with_columns((pl.col("end") - pl.col("start")).alias("weight"))
        .drop(["start", "end"])
        .collect()
        .to_pandas()
    )
    return df


df = polars_read()
print(df, file=sys.stderr)
coverage = weighted_median(df, "coverage", "weight")

min_coverage = get_min_coverage(coverage)
max_coverage = get_max_coverage(coverage)
mean = (df["coverage"] * df["weight"]).sum() / df["weight"].sum()
print(f"\nmean coverage: {mean}", file=sys.stderr)
print(f"median coverage: {coverage}\n", file=sys.stderr)

if coverage < 5:
    raise ValueError(
        f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. We recommend at least 10x coverage to use FIRE and require at least 5x."
        "If you are only examining data from a subset of chromosomes, consider using the keep_chromosomes parameter in config.yaml"
    )

open(arg.output_median, "w").write(str(round(coverage)) + "\n")
open(args.output_min, "w").write(str(round(min_coverage)) + "\n")
open(args.output_max, "w").write(str(round(max_coverage)) + "\n")
