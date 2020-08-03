singularity: "docker://continuumio/miniconda3"
include: "rules/common.smk"

rule targets:
    input:
        expand("{ref_fasta}.{suffix}", ref_fasta=ref_fasta, suffix=file_suffixes),
        "mafs/variants.maf",
        "qc/multiqc_report.html",
        "qc/depths.svg"

include: "rules/preprocessing.smk"
include: "rules/calling.smk"
include: "rules/pon.smk"
