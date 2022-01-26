# File: gather_depths.py
# Description: Aggregate exome depth statistics from mosdepth
# into a table for plotting

import pandas as pd

col_names = ['chr', 'start', 'end', 'depth']

def parse_summaries(files):
    depths = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        s = os.path.basename(f).split('.')
        patient = s[0]
        sample_type = s[1]
        depth = df.query("chrom == 'total_region'")
        depth['patient'] = patient
        depth['sample_type'] = sample_type
        depths.append(depth)
    return pd.concat(depths, ignore_index = True)
    
def main():
    files = snakemake.input
    df = parse_summaries(files)
    df = df[['patient', 'sample_type', 'chrom', 'length', 'bases', 'mean', 'min', 'max']]
    df.to_csv(snakemake.output[0], index=False, header=True)

if __name__ == '__main__':
    main()

