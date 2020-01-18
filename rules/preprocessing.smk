# File: preprocessing.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Holds rules to create an 'analysis ready bam' according
# to the GATK best practices guide. This bam can then be
# used for either germline or somatic variant calling


rule combine_fqs:
    input:
        unpack(get_fastq)
    output:
        temp("bams/{sample}.{type}.unaligned.bam")
    params:
        pl=get_platform
    shell:
        """
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} -O {output} \
            -SM {wildcards.sample}.{wildcards.type} \
            -RG {wildcards.sample}.{wildcards.type} \
            -PL {params.pl}
        """


rule bwa_index:
    input:
        ref_fasta
    output:
        [f"{ref_fasta}.{suffix}" for suffix in file_suffixes]
    shell:
        """
        bwa index {input}
        """

rule bwa:
    input:
        bam="bams/{sample}.{type}.unaligned.bam",
        ref=ref_fasta
    output:
        temp("bams/{sample}.{type}.aligned.bam")
    threads: 8
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
        unaligned="bams/{sample}.{type}.unaligned.bam",
        aligned="bams/{sample}.{type}.aligned.bam",
        ref=ref_fasta
    output:
        temp("bams/{sample}.{type}.merged.bam")
    shell:
        """
        gatk MergeBamAlignment -UNMAPPED {input.unaligned} -ALIGNED {input.aligned} \
            -R {input.ref} -O {output} 
        """

rule mark_duplicates:
    input:
        "bams/{sample}.{type}.merged.bam"
    output:
        bam=temp("bams/{sample}.{type}.markdups.bam"),
        md5=temp("bams/{sample}.{type}.markdups.bam.md5"),
        metrics="qc/gatk/{sample}_{type}_dup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} \
            --CREATE_MD5_FILE true --ASSUME_SORT_ORDER "queryname"
        """

rule sort_fix_tags:
    input:
        bam="bams/{sample}.{type}.markdups.bam",
        ref=ref_fasta
    output:
        bam=temp("bams/{sample}.{type}.sorted.bam"),
        bai=temp("bams/{sample}.{type}.sorted.bai"),
        md5=temp("bams/{sample}.{type}.sorted.bam.md5")
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
        bam="bams/{sample}.{type}.sorted.bam",
        known=known_sites,
        ref=ref_fasta
    output:
        bam="bams/{sample}.{type}.bam",
        bai="bams/{sample}.{type}.bai",
        md5="bams/{sample}.{type}.bam.md5",
        recal="qc/{sample}.{type}.recal_data.table"
    params:
        ks=['--known-sites ' + s for s in known_sites]
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} -O {output.recal} \
            {params.ks}
        gatk ApplyBQSR -I {input.bam} -R {input.ref} -O {output.bam} -bqsr {output.recal} \
            --static-quantized-quals 10 --static-quantized-quals 20 \
            --static-quantized-quals 30 --add-output-sam-program-record \
            --create-output-bam-md5 --use-original-qualities
        """

rule exome_cov:
    input:
        bam="bams/{sample}.{type}.bam",
        exons=capture_bed
    output:
        "qc/{sample}_{type}.mosdepth.region.dist.txt"
    threads: 4
    shell:
        """
        mosdepth --by {input.exons} -t {threads} qc/{wildcards.sample}_{wildcards.type} \
            {input.bam}
        """

rule stats:
    input:
        "bams/{sample}.{type}.bam"
    output:
        "qc/{sample}.{type}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule fastqc:
    input:
        "bams/{sample}.{type}.bam"
    output:
        html="qc/fastqc/{sample}_{type}.html",
        zip="qc/fastqc/{sample}_{type}_fastqc.zip"
    wrapper:
        "0.45.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_{type}_fastqc.zip", sample=samples, type=types),
        expand("qc/{sample}.{type}.recal_data.table", sample=samples, type=types),
        expand("qc/{sample}_{type}.mosdepth.region.dist.txt", sample=samples, type=types),
        expand("qc/{sample}.{type}.flagstat", sample=samples, type=types)
    output:
        "qc/multiqc_report.html"
    shell:
        """
        multiqc {input} -o qc/
        """
