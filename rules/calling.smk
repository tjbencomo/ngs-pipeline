# File: calling.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Contains rules to call somatic variants
# using GATK


rule mutect2:
    input:
        # unpack(get_call_pair),
        unpack(get_mutect2_input),
        ref=ref_fasta,
        germ_res=germline_resource
    output:
        vcf="vcfs/{patient}.unfiltered.vcf",
        idx="vcfs/{patient}.unfiltered.vcf.idx",
        stats="vcfs/{patient}.unfiltered.vcf.stats"
    params:
        tumor="{patient}.tumor",
        normal="{patient}.normal",
        pon=lambda wildcards, input: "--panel-of-normals " + pon_vcf,
        extra=""
    threads: 4
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.normal} -I {input.tumor} \
            -normal {params.normal} -tumor {params.tumor} -O {output.vcf} \
            --germline-resource {input.germ_res} \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            {params.pon} {params.extra}
        """

rule pileup_summaries:
    input:
        bam="bams/{patient}.tumor.bam",
        germ_res=contamination_resource
    output:
        "qc/{patient}_pileupsummaries.table"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} -V {input.germ_res} \
            -L {input.germ_res} -O {output}
        """

rule calculate_contamination:
    input:
        "qc/{patient}_pileupsummaries.table"
    output:
        "qc/{patient}_contamination.table"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk CalculateContamination -I {input} -O {output}
        """

rule filter_calls:
    input:
        vcf="vcfs/{patient}.unfiltered.vcf",
        ref=ref_fasta,
        contamination="qc/{patient}_contamination.table"
    output:
        vcf=temp("vcfs/{patient}.filtered.vcf.gz"),
        idx=temp("vcfs/{patient}.filtered.vcf.gz.tbi"),
        intermediate=temp("vcfs/{patient}.unselected.vcf"),
        inter_idx=temp("vcfs/{patient}.unselected.vcf.filteringStats.tsv"),
        inter_stats=temp("vcfs/{patient}.unselected.vcf.idx")
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
        vcf="vcfs/{patient}.filtered.vcf.gz",
        fasta=vep_fasta
    output:
        vcf="vcfs/{patient}.vcf",
        stats="vcfs/{patient}.vcf_summary.html"
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
        vcf="vcfs/{patient}.vcf",
        fasta=vep_fasta,
        vep_dir=vep_dir
    output:
        "mafs/{patient}.maf"
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
            --tumor-id {wildcards.patient}.tumor \
            --normal-id {wildcards.patient}.normal \
            --ref-fasta {input.fasta} --vep-data {input.vep_dir} \
            --ncbi-build {params.assembly} \
            --filter-vcf 0 --vep-path $vep_path \
            --maf-center {params.center}
        """
rule concat_mafs:
    input:
        expand("mafs/{patient}.maf", patient=patients)
    output:
        "mafs/variants.maf"
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/combine_mafs.py"
