include: "rules/common.smk"

rule targets:
    input:
        expand("bams/{sample}.{type}.bam", sample=samples, type=types),
        "qc/multiqc_report.html"

include: "rules/preprocessing.smk"
