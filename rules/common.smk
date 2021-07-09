# File: common.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Initialization file that loads config information and
# stores functions to aide with wildcard creation

import sys
import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.9.1")

configfile: "config.yaml"
validate(config, schema = "../schemas/config.schema.yaml")
patients = pd.read_csv(config['patients'])['patient']
units = pd.read_csv(config['units'], dtype=str).set_index(["patient", "sample", "readgroup"], drop=False)
validate(units, schema = "../schemas/units.schema.yaml")
units = units.sort_index()
seqtype= config['sequencing_type']

ref_dir = config['ref_dir']
ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
ref_dict = ref_fasta.replace('.fasta', '.dict') #used for BedToIntervals
known_sites = config['known_sites'].replace(' ', '').split(',')
known_sites = [os.path.join(ref_dir, s) for s in known_sites]
capture_bed = config['capture_targets']
germline_resource = config['germline_resource']
contamination_resource = config['contamination_resource']

num_workers = config['num_workers']

stringent_filtering = config['stringent_filtering']

tmp_dir = config['tmp_dir']
if tmp_dir == 'None':
    tmp_dir = 'null'

vep_dir = config['vep_dir']
assembly = config['assembly_version']
center = config['center_name']
alternate_isoforms = config['alternate_isoforms']
if config['alternate_isoforms'] == 'None':
    alternate_isoforms = ''
    isoforms_param = ''
else:
    isoforms_param = '--custom-enst ' + alternate_isoforms

tumor_only = config['tumor_only']

use_pon = config['use_pon']
if use_pon is False:
    pon_vcf = 'null'
else :
    if config['pon_vcf'] == 'None':
        build_pon = True 
        pon_vcf = "pon/pon.vcf.gz"
    else:
        build_pon = False
        pon_vcf = config['pon_vcf']

regions_bed = config['genomic_regions']
regions_gatk = os.path.basename(regions_bed).replace('.bed', '.interval_list')
regions_gatk = os.path.join('interval-files', regions_gatk)

if tumor_only:
    sample_types = ['tumor']
else:
    sample_types = ['normal', 'tumor']

file_suffixes= ['amb', 'ann', 'bwt', 'pac', 'sa']

# Load container info
gatk_env = config['gatk_container']
multiqc_env = config['multiqc_container']
mosdepth_env = config['mosdepth_container']
fastqc_env = config['fastqc_container']
vep_env = config['vep_container']
eda_env = config['eda_container']

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
    return units.loc[(wildcards.patient, wildcards.sample_type, wildcards.readgroup), 'platform']

def get_mutect2_input(wildcards):
    files = {}
    files['tumor'] = f"bams/{wildcards.patient}.tumor.bam"
    if use_pon:
        files['pon'] = pon_vcf
    if not tumor_only:
        files['normal'] = f"bams/{wildcards.patient}.normal.bam"
    return files
    # if use_pon:
    #     return {
    #         'normal' : f"bams/{wildcards.patient}.normal.bam",
    #         'tumor' : f"bams/{wildcards.patient}.tumor.bam",
    #         'pon' : pon_vcf
    #     }
    # else:
    #     return {
    #             'normal' : f"bams/{wildcards.patient}.normal.bam",
    #             'tumor' : f"bams/{wildcards.patient}.tumor.bam",
    #     }


def get_call_pair(wildcards):
    return {'normal' : f"bams/{wildcards.patient}.normal.bam",
            'tumor' : f"bams/{wildcards.patient}.tumor.bam"}

def get_vcf2maf_input(wildcards):
    if alternate_isoforms == '':
        return {
            'vcf' : f"vcfs/{wildcards.patient}.vcf",
            'fasta' : ref_fasta,
            'vep_dir' : vep_dir
        }
    else:
        return {
            'vcf' : f"vcfs/{wildcards.patient}.vcf",
            'fasta' : ref_fasta,
            'vep_dir' : vep_dir,
            'alt_isoforms' : alternate_isoforms
        }

def get_contamination_input(wildcards):
    out = {}
    out['tumor'] = f'qc/{wildcards.patient}_tumor_pileupsummaries.table'
    if not tumor_only:
        out['normal'] = f'qc/{wildcards.patient}_normal_pileupsummaries.table'
    return out

def get_coverage_input(wildcards):
    files = {}
    files['bam'] = f"bams/{wildcards.patient}.{wildcards.sample_type}.bam"
    if seqtype == "WES":
        files['regions'] = regions_bed
    return files

def isWGS(wildcards):
    seqtype = units.loc[(wildcards.patient, wildcards.sample_type), 'seqtype'][0]
    return seqtype == "WGS"


def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = '0' * num_zeros + str(i)
        ints.append(interval)
    return ints
   
def get_interval_files():
    ints = get_intervals()
    files = [i + '-scattered.interval_list' for i in ints]
    files = [os.path.join("interval-files", f) for f in files]
    return files

def get_orientationbias_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.patient}.{i}.f1r2.tar.gz" for i in intervals]
    return files

def get_mergevcfs_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.patient}.{i}.unfiltered.vcf" for i in intervals]
    return files

def get_mergestats_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.patient}.{i}.unfiltered.vcf.stats" for i in intervals]
    return files

interval_files = get_interval_files()
