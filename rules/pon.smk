# File: pon.smk
# Author: Tomas Bencomo
# Description: Workflow to create Panel of Normals for Mutect2

rule mutect2_pon:
    input:
        bam="bams/{patient}.normal.bam",
        ref=ref_fasta
    output:
        "pon/{patient}.pon.vcf.gz"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk Mutect2 -I {input.bam} -R {input.ref} -O {output} \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        -max-mnp-distance 0
        """
rule gather_variants:
    input:
        vcfs=expand("pon/{patient}.pon.vcf.gz", patient=patients),
        ref=ref_fasta,
        intervals=genome_intervals
    output:
        directory("pon/pon_db")
    params:
        vcfs=lambda wildcards, input: " -V ".join(input.vcfs)
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk GenomicsDBImport -R {input.ref} --genomicsdb-workspace-path {output} \
            -V {params.vcfs} -L {input.intervals}
        """

rule create_pon:
    input:
        var="pon/pon_db",
        ref=ref_fasta
    output:
        "pon/pon.vcf.gz"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk CreateSomaticPanelOfNormals -R {input.ref} -V gendb://{input.var} -O {output}
        """

