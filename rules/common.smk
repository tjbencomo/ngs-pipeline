# File: common.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Initialization file that loads config information and
# stores functions to aide with wildcard creation

import sys
import pandas as pd

configfile: "config.yaml"
samples = pd.read_csv(config['samples'])['sample']
units = pd.read_csv(config['units'], dtype=str).set_index("sample", drop=False)
ref_dir = config['ref_dir']
ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
known_sites = config['known_sites'].split(',')
known_sites = [os.path.join(ref_dir, s) for s in known_sites]
capture_bed = config['exome_targets']

wildcard_constraints:
    sample="|".join(samples)

def get_fastq(sample):
    return {'r1' : units.loc[sample, 'fq1'], 'r2' : units.loc[sample, 'fq2']}

def get_platform(sample):
    return units.loc[sample, 'platform'][0]


