# File: calling.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Contains rules to call somatic variants
# using GATK


rule mutect2:
    input:
        get_call_pair,
        ref=ref_fasta,
        germ_res=germline_resource
    output:
        vcf="vcfs/{sample}.unfiltered.vcf",
        bam="{bams/mutect2/{sample}_m2.bam}"
    params:
        tumor="{sample}.tumor}",
        normal="{sample}.normal",
        extra=""
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.normal} -I {input.tumor} \
            -normal {params.normal} -tumor {params.tumor} -O {output.vcf} \
            -bamout {output.bam}
        """



