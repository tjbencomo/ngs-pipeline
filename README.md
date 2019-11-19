# ngs-pipelines
Pipelines for processing next generation sequencing (NGS) data.

## Usage
Clone `ngs-pipelines` into the directory where your data is stored. After
specifying the `samples.csv` and `units.csv` run the desired pipeline with
```
conda activate snakemake
snakemake [pipeline]
```
`[pipeline]` can be either `preprocess`, `variant-call`, `rnaseq`.

## Installation
