singularity: "docker://continuumio/miniconda3"
include: "rules/common.smk"

rule targets:
    input:
        expand("bams/{sample}.{type}.bam", sample=samples, type=types),
        expand("vcfs/{sample}..vcf", sample=samples),
        expand("{ref_fasta}.{suffix}", ref_fasta=ref_fasta, suffix=file_suffixes),
        "qc/multiqc_report.html"

include: "rules/preprocessing.smk"
include: "rules/calling.smk"
