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
units = pd.read_csv(config['units'], dtype=str).set_index(["sample", "type"], drop=False)
ref_dir = config['ref_dir']
ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
known_sites = config['known_sites'].split(',')
known_sites = [os.path.join(ref_dir, s) for s in known_sites]
capture_bed = config['exome_targets']
types = ['control', 'tumor']

wildcard_constraints:
    sample="|".join(samples),
    type="|".join(types)

def get_fastq(wildcards):
    return {'r1' : units.loc[(wildcards.sample, wildcards.type), 'fq1'], 
            'r2' : units.loc[(wildcards.sample, wildcards.type), 'fq2']}

def get_platform(wildcards):
    return units.loc[(wildcards.sample, wildcards.type), 'platform'][0]

def get_call_pair(wildcards):
    return {'normal' : f"{wildcards.sample}.control.bam",
            'tumor' : f"{wildcards.sample}.tumor.bam"}

