"""
File: combine_mafs.py
Author: Tomas Bencomo
Email: tjbencomo@gmail.com
Description:
    This script concatenates a list of MAF formatted files
    into a single MAF file
"""

import os
import sys
import argparse
import pandas as pd

def get_args():
    """
    Not needed with snakemake handling input arguments, but useful if
    you're going to use this as a standalone script
    """
    parser = argparse.ArgumentParser(description="Combine MAF files")
    parser.add_argument('files', metavar='MAF', type=str, nargs='+',
                                help='Space separated list of MAF files')
    args = parser.parse_args()
    files = args.files
    for f in files:
        if not os.path.isfile(f):
            raise ValueError(f"File {f} is not a file or doesn't exist!")
    return args.files

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
