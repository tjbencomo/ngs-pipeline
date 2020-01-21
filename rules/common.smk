# File: common.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Initialization file that loads config information and
# stores functions to aide with wildcard creation

import sys
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.9.1")

configfile: "config.yaml"
validate(config, schema = "../schemas/config.schema.yaml")
samples = pd.read_csv(config['samples'])['sample']
units = pd.read_csv(config['units'], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema = "../schemas/units.schema.yaml")
ref_dir = config['ref_dir']
ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
known_sites = config['known_sites'].replace(' ', '').split(',')
known_sites = [os.path.join(ref_dir, s) for s in known_sites]
capture_bed = config['exome_targets']
germline_resource = config['germline_resource']
contamination_resource = config['contamination_resource']
vep_dir = config['vep_dir']
vep_fasta = config['vep_fasta']
assembly = config['assembly_version']
types = ['normal', 'tumor']
file_suffixes= ['amb', 'ann', 'bwt', 'pac', 'sa']

wildcard_constraints:
    sample="|".join(samples),
    type="|".join(types)

def get_fastq(wildcards):
    return {'r1' : units.loc[(wildcards.sample, wildcards.type), 'fq1'], 
            'r2' : units.loc[(wildcards.sample, wildcards.type), 'fq2']}

def get_platform(wildcards):
    return units.loc[(wildcards.sample, wildcards.type), 'platform'][0]

def get_call_pair(wildcards):
    return {'normal' : f"bams/{wildcards.sample}.normal.bam",
            'tumor' : f"bams/{wildcards.sample}.tumor.bam"}
