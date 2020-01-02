# ngs-pipeline
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
alias ngs='conda activate ngs-pipeline'
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
After finishing the setup and enabling the `conda` environment, inside the analysis directory with
`Snakefile` do a dry run to check for errors
```
snakemake -n
```
If you're processing lots of samples, the dry run may be slow from printing all the rules. Use
```
snakemake -n --quiet
```
Once you're ready to run the analysis type
```
snakemake
```
If your machine has multiple cores, you can exploit your computer's parallelism with
```
snakemake -j [cores]
```
This will run multiple rules simultaneously, speeding up the analysis.

After the analysis is complete don't forget to check `qc/multiqc_report.html` for
quality control results about the analysis.

### SLURM
There are two options to run the analysis on a SLURM cluster

1. Run snakemake process on node using `screen` and `nohup`
If you can reliably access a specific login node and there isn't a time limit
on login node processes, it's advisable to run the snakemake command that submits jobs
to the cluster on the login node with `screen` and `nohup`.
You can launch the analysis with cluster execution by typing
```
snakemake --cluster-config cluster.json -j 100 --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}  -c {cluster.ncpus}"
```
The benefit of this approach is `snakemake` temporary files are supported, so the analysis will
only create 4 files per sample: `{sample}.normal.bam`, `{sample}.tumor.bam`, `{sample}.vcf`, `{sample}_m2.bam`.
These files are the analysis ready bams, final filtered somatic calls, and the accompanying bam file 
for IGV viewing.

2. Submit all jobs to the queue immediately 
If your cluster has login node time limits or there are several login nodes, you may not be able to use
the first option. Time limits may kill the snakemake process prematurely. If there are several login nodes
and it varies which login node you're assigned to when you sign in, you may not be able to retrieve the process
if you launch it on one login node and then sign back in on another node.

To circumvent these problems, it's possible to submit all the jobs immediately to the queue,
returning control to the user. Unfortunately this approach doesn't support temporary files, so there will
be  intermediate files generated that you'll have to ignore or manually delete.


First give `parseJobID.sh` permission to run with
```
chmod +x parseJobID.sh
```

Then run snakemake with
```
snakemake --cluster-config cluster.json --cluster "sbatch $(./parseJobID.sh {dependencies}) -t {cluster.time} --mem {cluster.mem} -p {cluster.partition} -c {cluster.ncpus}" --jobs 100 --notemp --immediate-submit
```


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
