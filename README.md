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

1. Create a new Github repository using this workflow as a template using the `Use this template` button
at the top of this page. This will allow you to track any changes made to the analysis with `git`
2. Clone the repository to the machine where you want to perform data analysis
3. Create the `ngs-pipeline` 
environment with conda
```
conda env create -f environment.yml
```
This environment contains `snakemake` and the other executables (`samtools`, `gatk` etc) that you'll
need for data analysis.

4. Activate the environment with
```
conda activate ngs-pipeline
```
If you plan to use this environment frequently, its useful to create a bash alias for quick access.
Add this code to `.bashrc`
```
alias ngs='conda activate ngs-pipeline
```
After reloading `.bashrc` with `source .bashrc`, you can enable the environment by typing `ngs` at the console.

5. Configure `config.yml` to tell `ngs-pipeline` where to find important files for the workflow.

| Item                   | Description                                                                          |
|------------------------|--------------------------------------------------------------------------------------|
| samples                | CSV file with sample column. Each sample will have both a tumor and normal file      |
| units                  | CSV file the following columns: sample, type, platform, fq1, and fq2                 |
| ref_dir                | Filepath to directory where the reference directory is located                       |
| ref_fasta              | Name of the genome assembly file                                                     |
| known_sites            | Comma separated string with list of files with sites of known mutations              |
| exome_targets          | BED file with segments for coverage analysis                                         |
| germline_resource      | File with germline variants for Mutect2                                              |
| contamination_resource | File with biallelic germline variants for CalculateContamination                     |

`ref_fasta` and `known_sites` are expected to be in `ref_dir`. Only specify the filenames without paths.
In `ngs-pipeline`, each sample represents one patient. There should be normal and tumor sequencing data for each
sample. Each sample should have two rows in `units`, one normal row and one tumor row. Sequencing data must be
paired, so both `fq1` and `fq2` must be specified.


## Usage


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
