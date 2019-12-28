include: "rules/common.smk"

rule targets:
    input:
        expand("bams/{sample}.{type}.bam", sample=samples, type=types),
        expand("vcfs/{sample}.vcf", sample=samples),
        "qc/multiqc_report.html"

include: "rules/preprocessing.smk"
include: "rules/calling.smk"
