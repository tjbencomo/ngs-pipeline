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
import pandas as pd

def concat_mafs(files):
    mafs = []
    for f in files:
        df = pd.read_csv(f, sep = "\t", skiprows=1)
        mafs.append(df)
    all_mafs = pd.concat(mafs)
    all_mafs['patient'] = all_mafs['Tumor_Sample_Barcode'].str.replace(".tumor", "", regex=False)
    if all_mafs.shape[0] == 0:
        print("WARNING: There are no variants in your MAF file!")
    return all_mafs 

def main():
    files = snakemake.input
    maf = concat_mafs(files)
    if snakemake.params[0] is True:
        maf = maf.query('t_depth >= 15 & t_alt_count >= 5')
    maf.to_csv(snakemake.output[0], sep="\t", index=False, header=True)

if __name__ == '__main__':
    main()
