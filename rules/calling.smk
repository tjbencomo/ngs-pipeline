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
        vcf="vcfs/{sample}.unfiltered.vcf",
        idx="vcfs/{sample}.unfiltered.vcf.idx",
        stats="vcfs/{sample}.unfiltered.vcf.stats"
    params:
        tumor="{sample}.tumor",
        normal="{sample}.normal",
        extra=""
    conda:
        "../envs/gatk.yml"
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
    conda:
        "../envs/gatk.yml"
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
    conda:
        "../envs/gatk.yml"
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
        vcf=temp("vcfs/{sample}.filtered.vcf.gz"),
        idx=temp("vcfs/{sample}.filtered.vcf.gz.tbi"),
        intermediate=temp("vcfs/{sample}.unselected.vcf"),
        inter_idx=temp("vcfs/{sample}.unselected.vcf.filteringStats.tsv"),
        inter_stats=temp("vcfs/{sample}.unselected.vcf.idx")
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} -O {output.intermediate}
        gatk SelectVariants -V {output.intermediate} -R {input.ref} -O {output.vcf} \
            --exclude-filtered -OVI
        """

rule vep:
    input:
        vcf="vcfs/{sample}.filtered.vcf.gz",
        fasta=vep_fasta
    output:
        vcf="vcfs/{sample}.vcf",
        stats="vcfs/{sample}.vcf_summary.html"
    params:
        vep_dir = vep_dir,
        assembly=assembly
    conda:
        "../envs/annotation.yml"
    shell:
        """
        vep --cache --offline --hgvs --vcf --assembly {params.assembly} --dir {params.vep_dir}  \
            -i {input.vcf} -o {output.vcf} --fasta {input.fasta}
        """
rule vcf2maf:
    input:
        vcf="vcfs/{sample}.vcf",
        fasta=vep_fasta,
        vep_dir=vep_dir
    output:
        "mafs/{sample}.maf"
    conda:
        "../envs/annotation.yml"
    params:
        assembly=assembly,
        center=center
    shell:
        """
        vep_fp=`which vep`
        vep_path=$(dirname "$vep_fp")
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output} \
            --tumor-id {wildcards.sample}.tumor \
            --normal-id {wildcards.sample}.normal \
            --ref-fasta {input.fasta} --vep-data {input.vep_dir} \
            --ncbi-build {params.assembly} \
            --filter-vcf 0 --vep-path $vep_path \
            --center {params.center}
        """
