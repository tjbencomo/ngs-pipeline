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
        germ_res=germline_resource,
        interval="interval-files/{interval}-scattered.interval_list"
    output:
        vcf=temp("vcfs/{patient}.{interval}.unfiltered.vcf"),
        idx=temp("vcfs/{patient}.{interval}.unfiltered.vcf.idx"),
        stats=temp("vcfs/{patient}.{interval}.unfiltered.vcf.stats"),
        f1r2tar=temp("vcfs/{patient}.{interval}.f1r2.tar.gz")
    params:
        tumor="{patient}.tumor",
        normalname= ' ' if tumor_only else '-normal ' + "{patient}.normal",
        normal_input=lambda wildcards, input: ' ' if tumor_only else "-I " + input.normal,
        pon="--panel-of-normals " + pon_vcf,
        extra=mutect_flags
    singularity: gatk_env
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.tumor} \
            -tumor {params.tumor} -O {output.vcf} \
            --germline-resource {input.germ_res} \
            --f1r2-tar-gz {output.f1r2tar} \
            -L {input.interval} \
            -ip 100 \
            {params.normal_input} {params.normalname} \
            {params.pon} {params.extra}
        """

rule orientation_bias:
    input:
        get_orientationbias_input
    output:
        "vcfs/{patient}.read_orientation_model.tar.gz"
    params:
        i=lambda wildcards, input: ['-I ' + d for d in input]
    singularity: gatk_env
    shell:
        """
        gatk LearnReadOrientationModel {params.i} -O {output}
        """

rule pileup_summaries:
    input:
        bam="bams/{patient}.{sample_type}.bam",
        germ_res=contamination_resource
    output:
        "qc/{patient}_{sample_type}_pileupsummaries.table"
    singularity: gatk_env
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
    singularity: gatk_env
    shell:
        """
        gatk CalculateContamination -I {input.tumor}  \
            {params.matched} \
            -O {output}
        """

rule merge_vcfs:
    input:
        get_mergevcfs_input
    output:
        vcf="vcfs/{patient}.unfiltered.vcf"
    params:
        i=lambda wildcards, input: ['-I ' + vcf for vcf in input]
    singularity: gatk_env
    shell:
        """
        gatk MergeVcfs {params.i} -O {output.vcf}
        """

rule merge_stats:
    input:
        get_mergestats_input
    output:
        stats="vcfs/{patient}.unfiltered.vcf.stats"
    params:
        i=lambda wildcards, input: ['-stats ' + s for s in input]
    singularity: gatk_env
    shell:
        """
        gatk MergeMutectStats {params.i} -O {output.stats} 
        """

rule filter_calls:
    input:
        vcf="vcfs/{patient}.unfiltered.vcf",
        ref=ref_fasta,
        contamination="qc/{patient}_contamination.table",
        stats="vcfs/{patient}.unfiltered.vcf.stats",
        f1r2model="vcfs/{patient}.read_orientation_model.tar.gz"
    output:
        vcf="vcfs/filtered/{patient}.filtered.vcf",
        idx="vcfs/filtered/{patient}.filtered.vcf.idx",
        inter_stats="vcfs/filtered/{patient}.filtered.vcf.filteringStats.tsv",
    params:
        extra=filter_flags
    singularity: gatk_env
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} \
            --stats {input.stats} \
            -ob-priors {input.f1r2model} \
            -O {output.vcf} \
            {params.extra}
        """

rule select_calls:
    input:
        vcf="vcfs/filtered/{patient}.filtered.vcf",
        idx="vcfs/filtered/{patient}.filtered.vcf.idx",
        ref=ref_fasta
    output:
        vcf="vcfs/{patient}.vcf",
        idx="vcfs/{patient}.vcf.idx"
    singularity: gatk_env
    shell:
        """
        gatk SelectVariants -V {input.vcf} \
            -R {input.ref} \
            -O {output.vcf} \
            --exclude-filtered \
            -OVI
        """

rule vcf2maf:
    input:
        unpack(get_vcf2maf_input)
    output:
        vep_vcf="vcfs/{patient}.vep.vcf",
        maf="mafs/{patient}.maf"
    singularity: vep_env
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
            --ref-fasta {input.fasta} \
            --vep-data {input.vep_dir} \
            --ncbi-build {params.assembly} \
            --vep-path $vep_path \
            --maf-center {params.center} \
            {params.normalid} \
            {params.isoforms}
        """

rule concat_mafs:
    input:
        expand("mafs/{patient}.maf", patient=patients)
    output:
        "mafs/variants.maf"
    params:
        stringent_criteria = stringent_filtering
    singularity: eda_env
    script:
        "../scripts/combine_mafs.py"
