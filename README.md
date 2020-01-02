# ngs-pipelines
**CURRENTLY IN DEVELOPMENT**

NGS pipeline for calling somatic variants from whole exome sequencing (WES)
or whole genome sequencing (WGS) data

## Prerequisites
`conda` is required for `snakemake` and the bioinformatics executables.
If it's not already installed, go [here](https://www.anaconda.com/distribution/) 
to download and install Anaconda.

Before you start setting up the pipeline, make sure you have your reference genome assembled.
The GRCh38 (hg38) genome is available on the Broad's
GATK [website](https://software.broadinstitute.org/gatk/download/bundle). 

After downloading the assembly files, run `bwa index` on the FASTA file and save the
index files to the same directory where the assembly files are stored. `bwa mem` won't
run properly without these index files.

## Setup
First clone this repository and then create the `ngs-pipeline` 
environment with conda
```
conda env create -f environment.yml
```
This environment contains `snakemake` and all the bioinformatics tools (`samtools`, `gatk` etc)
needed for the workflow. You'll need to activate the workflow with `conda activate ngs-pipeline`
before you can use `snakemake`

## Usage
Clone `ngs-pipelines` into the directory where your data is stored. After
specifying the `samples.csv` and `units.csv` run the desired pipeline with
```
conda activate snakemake
snakemake [pipeline]
```

## Slurm Parallelism
https://hpc-carpentry.github.io/hpc-python/17-cluster/

## Citations
This pipeline is based on `dna-seq-gatk-variant-calling` by 
[Johannes Köster](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling).
If you use this pipeline, make sure to cite references for all of the tools used in the workflow:
```
snakemake
gatk
samtools
mosdepth
fastqc
multiqc
vep
```
Citations to be added...
