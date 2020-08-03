# File: calling.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Contains rules to call somatic variants
# using GATK


rule mutect2:
    input:
        unpack(get_mutect2_input),
        ref=ref_fasta,
        germ_res=germline_resource
    output:
        vcf="vcfs/{patient}.unfiltered.vcf",
        idx="vcfs/{patient}.unfiltered.vcf.idx",
        stats="vcfs/{patient}.unfiltered.vcf.stats",
        f1r2tar="vcfs/{patient}.f1r2.tar.gz"
    params:
        tumor="{patient}.tumor",
        normalname= ' ' if tumor_only else '-normal ' + "{patient}.normal",
        normal_input=lambda wildcards, input: ' ' if tumor_only else "-I " + input.normal,
        pon="--panel-of-normals " + pon_vcf,
        extra=""
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.tumor} \
            -tumor {params.tumor} -O {output.vcf} \
            --germline-resource {input.germ_res} \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --f1r2-tar-gz {output.f1r2tar} \
            {params.normal_input} {params.normalname} \
            {params.pon} {params.extra}
        """
        # """
        # gatk Mutect2 -R {input.ref} -I {input.normal} -I {input.tumor} \
        #     -normal {params.normal} -tumor {params.tumor} -O {output.vcf} \
        #     --germline-resource {input.germ_res} \
        #     --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        #     {params.pon} {params.extra}
        # """

rule orientation_bias:
    input:
        "vcfs/{patient}.f1r2.tar.gz"
    output:
        "vcfs/{patient}.read_orientation_model.tar.gz"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk LearnReadOrientationModel -I {input} -O {output}
        """

rule pileup_summaries:
    input:
        bam="bams/{patient}.{sample_type}.bam",
        germ_res=contamination_resource
    output:
        "qc/{patient}_{sample_type}_pileupsummaries.table"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} -V {input.germ_res} \
            -L {input.germ_res} -O {output}
        """

rule calculate_contamination:
    input:
        unpack(get_contamination_input)
    output:
        "qc/{patient}_contamination.table"
    params:
        matched=lambda wildcards, input:'' if tumor_only else '-matched ' + input.normal
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk CalculateContamination -I {input}  \
        {params.matched} \
        -O {output}
        """

rule filter_calls:
    input:
        vcf="vcfs/{patient}.unfiltered.vcf",
        ref=ref_fasta,
        contamination="qc/{patient}_contamination.table",
        stats="vcfs/{patient}.unfiltered.vcf.stats",
        f1r2model="vcfs/{patient}.read_orientation_model.tar.gz"
    output:
        vcf="vcfs/{patient}.vcf",
        idx="vcfs/{patient}.vcf.idx",
        intermediate=temp("vcfs/{patient}.unselected.vcf"),
        inter_idx=temp("vcfs/{patient}.unselected.vcf.filteringStats.tsv"),
        inter_stats=temp("vcfs/{patient}.unselected.vcf.idx")
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} \
            --stats {input.stats} \
            --ob-priors {input.f1r2model} \
            -O {output.intermediate}
        gatk SelectVariants -V {output.intermediate} -R {input.ref} -O {output.vcf} \
            --exclude-filtered -OVI
        """

# rule vep:
#     input:
#         vcf="vcfs/{patient}.filtered.vcf.gz",
#         fasta=vep_fasta
#     output:
#         vcf="vcfs/{patient}.vcf",
#         stats="vcfs/{patient}.vcf_summary.html"
#     params:
#         vep_dir = vep_dir,
#         assembly=assembly
#     conda:
#         "../envs/annotation.yml"
#     shell:
#         """
#         vep --cache --offline --hgvs --vcf --assembly {params.assembly} --dir {params.vep_dir}  \
#             -i {input.vcf} -o {output.vcf} --fasta {input.fasta}
#         """

rule vcf2maf:
    input:
        unpack(get_vcf2maf_input)
        #vcf="vcfs/{patient}.vcf",
        #fasta=ref_fasta,
        #vep_dir=vep_dir,
        #alt_isoforms=alternate_isoforms
    output:
        vep_vcf="vcfs/{patient}.vep.vcf",
        maf="mafs/{patient}.maf"
    conda:
        "../envs/annotation.yml"
    params:
        assembly=assembly,
        center=center,
        isoforms=isoforms_param,
        normalid=lambda wildcards, input: '' if tumor_only else "--normal-id " + wildcards.patient + ".normal"
    shell:
        """
        vep_fp=`which vep`
        vep_path=$(dirname "$vep_fp")
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} \
            --tumor-id {wildcards.patient}.tumor \
            --ref-fasta {input.fasta} --vep-data {input.vep_dir} \
            --ncbi-build {params.assembly} \
            --filter-vcf 0 --vep-path $vep_path \
            --maf-center {params.center} \
            {params.normalid} \
            {params.isoforms}
        """
        # """
        # vep_fp=`which vep`
        # vep_path=$(dirname "$vep_fp")
        # vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} \
        #     --tumor-id {wildcards.patient}.tumor \
        #     --normal-id {wildcards.patient}.normal \
        #     --ref-fasta {input.fasta} --vep-data {input.vep_dir} \
        #     --ncbi-build {params.assembly} \
        #     --filter-vcf 0 --vep-path $vep_path \
        #     --maf-center {params.center} \
        #     {params.isoforms}
        # """
rule concat_mafs:
    input:
        expand("mafs/{patient}.maf", patient=patients)
    output:
        "mafs/variants.maf"
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/combine_mafs.py"
