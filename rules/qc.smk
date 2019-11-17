rule fastqc:
    input:
        "data/samples/{sample}.fastq"
    output:
        "qc/fastqc/{sample}_fastqc.html",
        "qc/fastqc/{sample}_fastqc.zip"
    shell:
        """
        ml biology fastqc
        fastqc -o qc/fastqc/ {input}
        """

rule multiqc:
    input:
        fq=expand("qc/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "qc/multiqc_report.html"
    shell:
        """
        ml biology py-multiqc
        multiqc -o qc/ qc/fastqc
        """
