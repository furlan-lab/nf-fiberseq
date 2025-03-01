#!/usr/bin/env python
import re
import logging
import sys
import argparse
import os
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Compute coverage statistics")
    parser.add_argument("--chromosomes",required=True, help="Txt file with chromosomes list")
    return parser.parse_args()

args = get_args()

def read_chromosomes():
    with open(args.chromosomes, "r") as f:
        chromosomes = [line.strip() for line in f if line.strip()]  # Remove empty lines
    return chromosomes

chroms = read_chromosomes()
print(chroms)