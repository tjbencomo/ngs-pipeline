# File: calling.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Contains rules to call somatic variants
# using GATK


rule mutect2:
    input:
        unpack(get_call_pair),
        ref=ref_fasta,
        germ_res=germline_resource
    output:
        vcf=temp("vcfs/{sample}.unfiltered.vcf")
    params:
        tumor="{sample}.tumor",
        normal="{sample}.normal",
        extra=""
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.normal} -I {input.tumor} \
            -normal {params.normal} -tumor {params.tumor} -O {output.vcf} \
            --germline-resource {input.germ_res} \
            {params.extra}
        """

rule pileup_summaries:
    input:
        bam="bams/{sample}.tumor.bam",
        germ_res=contamination_resource
    output:
        "qc/{sample}_pileupsummaries.table"
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} -V {input.germ_res} \
            -L {input.germ_res} -O {output}
        """

rule calculate_contamination:
    input:
        "qc/{sample}_pileupsummaries.table"
    output:
        "qc/{sample}_contamination.table"
    shell:
        """
        gatk CalculateContamination -I {input} -O {output}
        """

rule filter_calls:
    input:
        vcf="vcfs/{sample}.unfiltered.vcf",
        ref=ref_fasta,
        contamination="qc/{sample}_contamination.table"
    output:
        vcf="vcfs/{sample}.vcf.gz",
        idx="vcfs/{sample}.vcf.gz.tbi",
        intermediate=temp("vcfs/{sample}.unselected.vcf")
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} -O {output.intermediate}
        gatk SelectVariants -V {output.intermediate} -R {input.ref} -O {output.vcf} \
            --exclude-filtered -OVI
        """




