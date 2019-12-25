# ngs-pipelines
Pipelines for processing next generation sequencing (NGS) data.

## Prerequisites
`conda` is required for `snakemake` and the bioinformatics executables.
If it's not already installed, go [here](https://www.anaconda.com/distribution/) 
to download and install Anaconda.

## Setup
First clone this repository and then create the `ngs-pipeline` 
environment with conda
```
conda env create -f environment.yml
```

## Usage
Clone `ngs-pipelines` into the directory where your data is stored. After
specifying the `samples.csv` and `units.csv` run the desired pipeline with
```
conda activate snakemake
snakemake [pipeline]
```
`[pipeline]` can be either `preprocess`, `variant-call`, `rnaseq`.

## Installation
