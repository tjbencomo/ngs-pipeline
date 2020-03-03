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
        temp("bams/{patient}.{sample_type}.{readgroup}.unaligned.bam")
    params:
        pl=get_platform
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} -O {output} \
            -SM {wildcards.patient}.{wildcards.sample_type} \
            -RG {wildcards.patient}.{wildcards.sample_type}.{wildcards.readgroup} \
            -PL {params.pl}
        """


rule bwa_index:
    input:
        ref_fasta
    output:
        [f"{ref_fasta}.{suffix}" for suffix in file_suffixes]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        bwa index {input}
        """

rule bwa:
    input:
        bam="bams/{patient}.{sample_type}.{readgroup}.unaligned.bam",
        ref=ref_fasta
    output:
        temp("bams/{patient}.{sample_type}.{readgroup}.aligned.bam")
    threads: 8
    conda:
        "../envs/gatk.yml"
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
        unaligned="bams/{patient}.{sample_type}.{readgroup}.unaligned.bam",
        aligned="bams/{patient}.{sample_type}.{readgroup}.aligned.bam",
        ref=ref_fasta
    output:
        temp("bams/{patient}.{sample_type}.{readgroup}.merged.bam")
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeBamAlignment -UNMAPPED {input.unaligned} -ALIGNED {input.aligned} \
            -R {input.ref} -O {output} 
        """

rule mark_duplicates:
    input:
        get_dedup_input
        # "bams/{patient}.{sample_type}.merged.bam"
    output:
        bam=temp("bams/{patient}.{sample_type}.markdups.bam"),
        md5=temp("bams/{patient}.{sample_type}.markdups.bam.md5"),
        metrics="qc/gatk/{patient}_{sample_type}_dup_metrics.txt"
    params:
        input=lambda wildcards, input: " -I  ".join(input)
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MarkDuplicates -I {params.input} -O {output.bam} -M {output.metrics} \
            --CREATE_MD5_FILE true --ASSUME_SORT_ORDER "queryname"
        """

rule sort_fix_tags:
    input:
        bam="bams/{patient}.{sample_type}.markdups.bam",
        ref=ref_fasta
    output:
        bam=temp("bams/{patient}.{sample_type}.sorted.bam"),
        bai=temp("bams/{patient}.{sample_type}.sorted.bai"),
        md5=temp("bams/{patient}.{sample_type}.sorted.bam.md5")
    conda:
        "../envs/gatk.yml"
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
        bam="bams/{patient}.{sample_type}.sorted.bam",
        known=known_sites,
        ref=ref_fasta
    output:
        bam="bams/{patient}.{sample_type}.bam",
        bai="bams/{patient}.{sample_type}.bai",
        md5="bams/{patient}.{sample_type}.bam.md5",
        recal="qc/{patient}.{sample_type}.recal_data.table"
    params:
        ks=['--known-sites ' + s for s in known_sites]
    conda:
        "../envs/gatk.yml"
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
        bam="bams/{patient}.{sample_type}.bam",
        exons=capture_bed
    output:
        "qc/{patient}_{sample_type}.mosdepth.region.dist.txt"
    threads: 4
    conda:
        "../envs/qc.yml"
    shell:
        """
        mosdepth --by {input.exons} -t {threads} qc/{wildcards.patient}_{wildcards.sample_type} \
            {input.bam}
        """

rule stats:
    input:
        "bams/{patient}.{sample_type}.bam"
    output:
        "qc/{patient}.{sample_type}.flagstat"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule fastqc:
    input:
        "bams/{patient}.{sample_type}.bam"
    output:
        html="qc/fastqc/{patient}_{sample_type}.html",
        zip="qc/fastqc/{patient}_{sample_type}_fastqc.zip"
    conda:
        "../envs/qc.yml"
    wrapper:
        "0.45.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{patient}_{sample_type}_fastqc.zip", patient=patients, sample_type=sample_types),
        expand("qc/{patient}_{sample_type}.mosdepth.region.dist.txt", patient=patients, sample_type=sample_types),
        expand("qc/{patient}.{sample_type}.flagstat", patient=patients, sample_type=sample_types)
    output:
        "qc/multiqc_report.html"
    conda:
        "../envs/qc.yml"
    shell:
        """
        multiqc {input} -o qc/
        """
