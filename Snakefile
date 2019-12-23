import os
import sys
import pandas as pd

configfile: "config.yaml"
samples = pd.read_csv(config['samples'])['sample']
units = pd.read_csv(config['units'], dtype=str).set_index("sample", drop=False)
ref_dir = config['ref_dir']
ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
known_sites = config['known_sites'].split(',')
known_sites = [os.path.join(ref_dir, s) for s in known_sites]

def get_fastq(sample):
    return {'r1' : units.loc[sample, 'fq1'], 'r2' : units.loc[sample, 'fq2']}

rule targets:
    input:
        expand("bams/{sample}.bam", sample=samples),
        "qc/multiqc_report.html"

rule combine_fqs:
    input:
        unpack(get_fastq)
    output:
        "bams/{sample}.unaligned.bam"
    shell:
        """
        ml biology gatk
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} -O {output} \
            -SM {wildcards.sample} -RG {wildcards.sample}
        """

rule bwa:
    input:
        bam="bams/{sample}.unaligned.bam",
        ref=ref_fasta
    output:
        "bams/{sample}.aligned.bam"
    threads: 4
    shell:
        """
        ml biology bwa samtools gatk
        gatk SamToFastq -I {input.bam} -F /dev/stdout -INTER true -NON_PF true \
        | \
        bwa mem -K 100000000 -p -v 3 -t {threads} -Y \
            {input.ref} /dev/stdin - 2> >(tee {log} >&2) \
        | \
        samtools view -1 - > {output}
        """

rule merge_bams:
    input:
        unaligned="bams/{sample}.unaligned.bam",
        aligned="bams/{sample}.aligned.bam",
        ref=ref_fasta
    output:
        "bams/{sample}.merged.bam"
    shell:
        """
        ml biology gatk
        gatk MergeBamAlignment -UNMAPPED {input.unaligned} -ALIGNED {input.aligned} \
            -R {input.ref} -O {output} 
        """

rule mark_duplicates:
    input:
        "bams/{sample}.merged.bam"
    output:
        bam="bams/{sample}.markdups.bam",
        md5="bams/{sample}.markdups.bam.md5",
        metrics="qc/gatk/{sample}_dup_metrics.txt"
    shell:
        """
        ml biology gatk
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} \
            --CREATE_MD5_FILE true --ASSUME_SORT_ORDER "queryname"
        """

rule sort_fix_tags:
    input:
        bam="bams/{sample}.markdups.bam",
        ref=ref_fasta
    output:
        bam="bams/{sample}.sorted.bam",
        bai="bams/{sample}.sorted.bai"
    shell:
        """
        ml biology gatk
        gatk SortSam -I {input.bam} -O /dev/stdout --SORT_ORDER "coordinate" \
            --CREATE_INDEX false --CREATE_MD5_FILE false \
        | \
        gatk SetNmMdAndUqTags -I /dev/stdin -O {output.bam} -R {input.ref} \
            --CREATE_INDEX true --CREATE_MD5_FILE true
        """

rule compute_bqsr:
    input:
        bam="bams/{sample}.sorted.bam",
        known=known_sites,
        ref=ref_fasta
    output:
        "qc/{sample}.recal_data.table"
    shell:
        """
        ml biology gatk
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} -O {output} \
            --known-sites {known} 
        """

rule apply_bqsr:
    input:
        bam="bams/{sample}.sorted.bam",
        ref=ref_fasta,
        bqsr="qc/{sample}.recal_data.table"
    output:
        "bams/{sample}.bam"
    shell:
        """
        ml biology gatk
        gatk ApplyBQSR -I {input.bam} -R {input.ref} -O {output} -bqsr {input.bqsr} \
            --static-quantized-quals 10 --static-quantized-quals 20 \
            --static-quantized-quals 30 --add-output-sam-program-record \
            --create-output-bam-md5 --use-original-qualities
        """

rule fastqc:
    input:
        "bams/{sample}.markdups.bam"
    output:
        "qc/fastqc/{sample}_fastqc.html",
        "qc/fastqc/{sample}_fastqc.zip"
    shell:
        """
        ml biology fastqc
        fastqc {input} -o qc/fastqc
        """

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report.html"
    shell:
        """
        ml biology py-multiqc
        multiqc {input} -o qc/
        """
