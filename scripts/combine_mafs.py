"""
File: combine_mafs.py
Author: Tomas Bencomo
Email: tjbencomo@gmail.com
Description:
    This script concatenates a list of MAF formatted files
    into a single MAF file. Designed to aggregate variants
    from many samples into a single analysis-wide MAF file.
"""

import os
import sys
import argparse
import pandas as pd

def concat_mafs(files):
    mafs = []
    for f in files:
        df = pd.read_csv(f, sep = "\t", skiprows=1)
        mafs.append(df)
    all_mafs = pd.concat(mafs)
    all_mafs['patient'] = all_mafs['Tumor_Sample_Barcode'].str.split(".", expand=True)[0]
    return all_mafs 

def main():
    files = snakemake.input
    maf = concat_mafs(files)
    maf.to_csv(snakemake.output[0], index=False)

if __name__ == '__main__':
    main()
