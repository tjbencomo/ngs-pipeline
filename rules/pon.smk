# File: pon.smk
# Author: Tomas Bencomo
# Description: Workflow to create Panel of Normals for Mutect2

rule mutect2_pon:
    input:
        bam="bams/{patient}.normal.bam",
        intervals=regions_gatk,
        ref=ref_fasta
    output:
        vcf="pon/{patient}.pon.vcf.gz",
        idx="pon/{patient}.pon.vcf.gz.tbi",
        stats="pon/{patient}.pon.vcf.gz.stats"
    singularity: gatk_env
    shell:
        """
        gatk Mutect2 \
            -I {input.bam} \
            -R {input.ref} \
            -O {output.vcf} \
            -L {input.intervals} \
            -ip 100 \
            --max-mnp-distance 0
        """
rule gather_variants:
    input:
        vcfs=expand("pon/{patient}.pon.vcf.gz", patient=patients),
        intervals=regions_gatk,
        ref=ref_fasta
    output:
        directory("pon/pon_db")
    params:
        vcfs=lambda wildcards, input: " -V ".join(input.vcfs)
    singularity: gatk_env
    shell:
        """
        gatk GenomicsDBImport \
            -R {input.ref} \
            --genomicsdb-workspace-path {output} \
            -V {params.vcfs} \
            -L {input.intervals} \
            -ip 100
        """

rule create_pon:
    input:
        var="pon/pon_db",
        ref=ref_fasta
    output:
        vcf="pon/pon.vcf.gz",
        idx="pon/pon.vcf.gz.tbi"
    singularity: gatk_env
    shell:
        """
        gatk CreateSomaticPanelOfNormals -R {input.ref} -V gendb://{input.var} -O {output.vcf}
        """

