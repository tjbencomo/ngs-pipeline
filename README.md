# ngs-pipeline
Somatic variant calling pipeline for whole exome and whole genome sequencing (WES & WGS). The pipeline uses Mutect2 to identify variants and mostly follows GATK Best Practices. SLURM execution functionality allows the workflow to run on Stanford's Sherlock computing cluster.

## Author
* Tomas Bencomo ([https://tjbencomo.github.io](https://tjbencomo.github.io))

## Getting Help
If you encounter problems using the pipeline, find a bug, or would like to request a new feature, 
please file an [issue](https://github.com/tjbencomo/ngs-pipeline/issues).

## Prerequisites
## References
You'll need a reference genome. The GRCh38 (hg38) genome is available on the Broad's
GATK [website](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

You'll also need the `cache` files for 
[Variant Annotation Predictor (VEP)](https://github.com/Ensembl/ensembl-vep).
Follow the tutorial 
[here](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) 
to download the data files. 
By default, the pipeline uses a docker container with VEP v104 and VCF2MAF for variant annotation. 
It is suggested you use v104 cache files. You can use your own version of VEP/VCF2MAF by
modifying the `vep_env` field in `config.yaml`.

## Input Files
The pipeline takes paired end (PE) FASTQ files as input. Samples can be normal-tumor pairs or
tumor-only individual sample. Samples can be divided into multiple read groups (see `units.csv`) or one pair
of FASTQ files per sample.

You'll also need a BED file with a list of target regions. For WES this is usually the capture kit regions. For WGS
this can be the entire genome or a whitelist file that excludes problematic regions. A WGS calling region file is available
in the GATK Resource Bundle (it will need to be converted from interval_list format to BED format).

**NOTE** If you have WES and WGS samples to analyze, create two separate instances of the workflow
and run the samples separately. 

## Software
Snakemake is required to run the pipeline. It is recommended users have Singularity installed to
take advantage of preconfigured Docker containers for full reproducibility. If you
don't want to use Singularity, you should download all required software (see workflow files) and
ensure they are in your path.

## Setup
1. Clone this repository or create a new repository using this workflow as a template
3. Edit `patients.csv` and `units.csv` with the details for your analysis.
See the `schemas/` directory for details about each file.
4. Configure `config.yml`. See `schemas/config.schema.yaml` for info about each required field. 
Multiplexed samples should be differentiated with the `readgroup` column.
Sequencing data must be paired, so both `fq1` and `fq2` are required.

## Usage
After finishing the setup, inside the repo's base  directory with
`Snakefile` do a dry run to check for errors
```
snakemake -n
```
Once you're ready to run the analysis type
```
snakemake --use-singularity -j [cores]
```
This will run the workflow locally. It is recommended you have at least 20GB of storage
available as the Singularity containers take up around 8GB and the output files can be very
large depending on the number of samples and sequencing depth.

The pipeline produces two key files: `mafs/variants.maf` and `qc/multiqc_report.html`.

`variants.maf` includes somatic variants from all samples that passed Mutect2 filtering.
They have been annotated with VEP and a single effect has been chosen by [vcf2maf](https://github.com/mskcc/vcf2maf)
using the Ensembl database. Ensembl uses its canonical isoforms for effect selection. 
Other isoforms can be specified with the `alternate_isoforms` field in `config.yaml`.
See the [cBioPortal override isoforms](https://github.com/mskcc/vcf2maf/blob/master/data/isoform_overrides_uniprot)
for file formatting.

`multiqc_report.html` includes quality metrics like coverage for the fully processed BAM files. 
Individual VCF files for each sample prior VEP annotation are found as `vcfs/{patient}.vcf`.
VEP annotated VCFs are found as `vcfs/{patient}.vep.vcf`. `qc/depths.svg` shows the sequencing depth distribution
for normal and tumor samples.

### Mutect2 Parallelism
Mutect2 can be scattered into many workers to process different regions of each sample in parallel via the `num_workers` setting in `config.yaml`.
This divides the `genomic_regions` file into many small subregions (the same number of subregions as the number of workers). 
Sets of intervals are split between files and individual intervals are not broken. When using the Broad's WGS calling regions interval file,
GATK does not seem to support more than 24 workers (GATK won't create more than 24 subregions). For WES data such as exome capture
kits, GATK will split the regions into more subsets (I have tested up to 50 workers but it can probably do more). 


### Cluster Execution
If you're using a compute cluster, you can take advantage of massively
parallel computation to speed up the analysis. Only SLURM clusters are
currently supported, but if you work with another cluster system (SGE etc)
Snakemake makes it relatively easy to add support for your cluster.

Follow these instructions to use SLURM execution
1. Modify `cluster.wes.json` if you are analyzing WES data or `cluster.wgs.json` for WGS data. Set the `account` field to your SLURM account and the `out` field to where the SLURM log files should be saved. The recommended format is `/path/slurm-{jobid}.out`. All directories in the path must exist before launch - Snakemake will not create any directories. 
2. Edit `run_pipeline.sh`. Specify the SLURM directives and set Snakemake to use `cluster.wes.json` or `cluster.wgs.json` depending on your needs.
3. When ready, launch the workflow
```
sbatch run_pipeline.sh [path to Snakemake directory]
```

If for some reason you can't leave the master `snakemake` process running, `snakemake`
offers the ability to launch all jobs using `--immediate-submit`. This
approach will submit all jobs to the queue immediately and finish the master `snakemake`
process. The downside of this method is that temporary files are not supported, so
the results directories will be very cluttered. 
To use `--immediate-submit` follow these steps
1. Give `parseJobID.sh` permission to run
```
chmod +x parseJobID.sh
```
2. Submit the `snakemake` jobs
```
snakemake --cluster-config cluster.json --cluster 'sbatch $(./parseJobID.sh {dependencies}) -t {cluster.time} --mem {cluster.mem} -p {cluster.partition} -c {cluster.ncpus} - o {cluster.out}' --jobs 100 --notemp --immediate-submit
```

## Deviations from Broad Pipeline
`ngs-pipeline` differs from the Broad's Somatic Variant pipeline in the following ways:
* By default `FilterMutectCalls` is run with additional flags. This can be changed in `config.yaml`. The default configuration requires all variants must have 1 read from each direction to be included in the final callset. 
* We provide the option for more stringent variant filtering criteria with the `stringent_filtering` setting in `config.yaml`. 
This is turned off by default


## Test Dataset
A small sample of `chr21` reads are supplied from the 
[Texas Cancer Research Biobank](http://txcrb.org/index.html) for
end-to-end pipeline tests. 
This data is made available as open access data with minimal privacy
restrictions. Please read the [Conditions of Use](http://txcrb.org/data.html)
before using the data.

## Tips
### FASTQ Formatting
The first line of each read (called the sequence identifier) should consist of a single string not separated
by any spaces. An example of a properly sequence identifier is
```
@MGILLUMINA4_74:7:1101:10000:100149
```
Sequence identifiers that consist of multiple space separated words will cause problems because `gatk` captures the entire
string and uses it as the read ID but `bwa` only parses the first word as the read group when it writes aligned reads
to BAM files. This sequence identifier would break the pipeline
```
@MGILLUMINA4_74:7:1101:10000:100149     RG:Z:A470018/1
```
Check the sequence identifiers if you encounter a `Aligned record iterator is behind the unmapped reads` error.

## Citations
This pipeline is based on `dna-seq-gatk-variant-calling` by 
[Johannes Köster](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling).
If you use `ngs-pipeline`, please use `citations.md` to cite the necessary software tools. 
