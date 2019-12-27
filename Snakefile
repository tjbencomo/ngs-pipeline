include: "common.smk"

rule targets:
    input:
        expand("bams/{sample}.bam", sample=samples),
        "qc/multiqc_report.html"

include: "preprocessing.smk"
