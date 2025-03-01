#!/usr/bin/env python

import os
import csv

# Change these for real data
# for this project - location of the hifi bam files are in: 
# /fh/fast/furlan_s/SR/ngs/pacbio/20231218_Mendoza_M/Fiberseq_Pool2.MM/r84046_20231208_203405_2_B01/hifi_reads/m84046_231209_111623_s2.hifi_reads.bc2013.bam
base_dir = "/fh/fast/furlan_s/SR/ngs/pacbio/20231218_Mendoza_M/"
sampleIdMap = {
    'bc2010': 'THP1',
    'bc2011': 'Kasumi',
    'bc2012': 'HL60',
    'bc2013': 'K562'
}
output_file = "/fh/fast/furlan_s/user/tho/nextflow/samplesheet.csv"

# Collect file paths and map to IDs
samples = []
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".bam") and "hifi_reads" in root:
            bam_path = os.path.join(root, file)
            for bc, sample_id in sampleIdMap.items():
                if bc in file:
                    samples.append({"id": sample_id, "reads_bam": bam_path})
                    break

# Write to CSV
with open(output_file, mode="w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["id", "reads_bam"])
    writer.writeheader()
    writer.writerows(samples)

print(f"Samplesheet generated: {output_file}")

