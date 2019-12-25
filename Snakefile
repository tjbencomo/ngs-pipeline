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
capture_bed = config['exome_targets']

wildcard_constraints:
    sample="|".join(samples)

def get_fastq(sample):
    return {'r1' : units.loc[sample, 'fq1'], 'r2' : units.loc[sample, 'fq2']}

def get_platform(sample):
    return units.loc[sample, 'platform'][0]

rule targets:
    input:
        expand("bams/{sample}.bam", sample=samples),
        "qc/multiqc_report.html"

rule combine_fqs:
    input:
        unpack(get_fastq)
    output:
        temp("bams/{sample}.unaligned.bam")
    params:
        pl=get_platform
    shell:
        """
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} -O {output} \
            -SM {wildcards.sample} -RG {wildcards.sample} -PL {params.pl}
        """

rule bwa:
    input:
        bam="bams/{sample}.unaligned.bam",
        ref=ref_fasta
    output:
        temp("bams/{sample}.aligned.bam")
    threads: 4
    shell:
        """
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
        temp("bams/{sample}.merged.bam")
    shell:
        """
        gatk MergeBamAlignment -UNMAPPED {input.unaligned} -ALIGNED {input.aligned} \
            -R {input.ref} -O {output} 
        """

rule mark_duplicates:
    input:
        "bams/{sample}.merged.bam"
    output:
        bam=temp("bams/{sample}.markdups.bam"),
        md5=temp("bams/{sample}.markdups.bam.md5"),
        metrics="qc/gatk/{sample}_dup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} \
            --CREATE_MD5_FILE true --ASSUME_SORT_ORDER "queryname"
        """

rule sort_fix_tags:
    input:
        bam="bams/{sample}.markdups.bam",
        ref=ref_fasta
    output:
        bam=temp("bams/{sample}.sorted.bam"),
        bai=temp("bams/{sample}.sorted.bai")
        md5=temp("bams/{sample}.sorted.bam.md5")
    shell:
        """
        gatk SortSam -I {input.bam} -O /dev/stdout --SORT_ORDER "coordinate" \
            --CREATE_INDEX false --CREATE_MD5_FILE false \
        | \
        gatk SetNmMdAndUqTags -I /dev/stdin -O {output.bam} -R {input.ref} \
            --CREATE_INDEX true --CREATE_MD5_FILE true
        """

rule bqsr:
    input:
        bam="bams/{sample}.sorted.bam",
        known=known_sites,
        ref=ref_fasta
    output:
        bam="bams/{sample}.bam",
        recal="qc/{sample}.recal_data.table"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} -O {output.recal} \
            --known-sites {input.known}
        gatk ApplyBQSR -I {input.bam} -R {input.ref} -O {output.bam} -bqsr {output.recal} \
            --static-quantized-quals 10 --static-quantized-quals 20 \
            --static-quantized-quals 30 --add-output-sam-program-record \
            --create-output-bam-md5 --use-original-qualities
        """

rule exome_cov:
    input:
        bam="bams/{sample}.bam",
        exons=capture_bed
    output:
        "qc/{sample}.mosdepth.region.dist.txt"
    threads: 4
    shell:
        """
        mosdepth --by {input.exons} -t {threads} qc/{wildcards.sample} {input.bam}
        """

rule stats:
    input:
        "bams/{sample}.bam"
    output:
        "qc/{sample}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule fastqc:
    input:
        "bams/{sample}.bam"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    wrapper:
        "0.45.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_fastqc.zip", sample=samples),
        expand("qc/{sample}.recal_data.table", sample=samples),
        expand("qc/{sample}.mosdepth.region.dist.txt", sample=samples),
        expand("qc/{sample}.flagstat", sample=samples)
    output:
        "qc/multiqc_report.html"
    shell:
        """
        multiqc {input} -o qc/
        """
