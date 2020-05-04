# File: count_depth.py
# Description: Aggregate exome depth statistics from mosdepth
# into a table for plotting

import pandas as pd

col_names = ['chr', 'start', 'end', 'depth']

def average_depth(files):
    patients = []
    samples = []
    depths = []
    for f in files:
        df = pd.read_csv(f, names = col_names, sep='\t')
        info = f.split('.')[0]
        patient = info.split('_')[0]
        sample_type = info.split('_')[1]
        patients.append(patient)
        samples.append(sample_type)
        depths.append(df['depth'].mean())
    return pd.DataFrame({'patient' : patients, 'type' : samples, 'mean_depth' : depths})

def main():
    files = snakemake.input
    df = average_depth(files)
    df.to_csv(snakemake.output[0], index=False, header=True)

if __name__ == '__main__':
    main()

