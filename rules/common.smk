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
patients = pd.read_csv(config['patients'])['patient']
units = pd.read_csv(config['units'], dtype=str).set_index(["patient", "sample", "readgroup"], drop=False)
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
center = config['center_name']
sample_types = ['normal', 'tumor']
file_suffixes= ['amb', 'ann', 'bwt', 'pac', 'sa']

wildcard_constraints:
    patient="|".join(patients),
    sample_type="|".join(sample_types)

def get_fastq(wildcards):
    return {'r1' : units.loc[(wildcards.patient, wildcards.sample_type, wildcards.readgroup), 'fq1'], 
            'r2' : units.loc[(wildcards.patient, wildcards.sample_type, wildcards.readgroup), 'fq2']}

def get_readgroups(wildcards):
    return units.loc[(wildcards.patient, wildcards.sample_type), 
                    'readgroup'].unique().tolist()

def get_dedup_input(wildcards):
    rgs = units.loc[(wildcards.patient, wildcards.sample_type), 'readgroup'].unique().tolist()
    bams = [f"bams/{wildcards.patient}.{wildcards.sample_type}.{rg}.merged.bam" for rg in rgs]
    return bams

def get_platform(wildcards):
    return units.loc[(wildcards.patient, wildcards.sample_type, wildcards.readgroup), 'platform'][0]

def get_call_pair(wildcards):
    return {'normal' : f"bams/{wildcards.patient}.normal.bam",
            'tumor' : f"bams/{wildcards.patient}.tumor.bam"}
