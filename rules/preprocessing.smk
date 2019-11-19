rule fastq_to_bam:
    input:
        "data/samples/{sample}.fastq"
    output:
        "results/{sample}_unaligned.bam"
    log:
        "logs/fastq_to_bam/{sample}.log"
    shell:
        """
        ml biology gatk
        gatk FastqToSam \
            --FASTQ={input} \
            --OUTPUT={output} \
            --PLATFORM="ILLUMINA" \
            --SAMPLE_NAME={wildcards.sample} \
            --READ_GROUP_NAME={wildcards.sample} 2> {log}
        """

rule bwa:
    input:
        bam="results/{sample}_unaligned.bam",
        ref="data/genome.fa"
    output:
        aligned="results/{sample}_aligned.bam",
    log:
        "logs/bwa/{sample}.log"
    threads: 4
    shell:
        """
        ml biology bwa samtools gatk
        gatk SamToFastq \
            --INPUT={input.bam} \
            --FASTQ=/dev/stdout \
            --INTERLEAVE=true \
            --INCLUDE_NON_PF_READS=true \
            | \
        bwa mem -K 100000000 -p -v 3 -t {threads} -Y \
            {input.ref} /dev/stdin - 2> >(tee {log} >&2) \
            | \
        samtools view -1 - > {output.aligned}
        """

rule merge_bam_alignment:
    input:
        aligned="results/{sample}_unaligned.bam",
        unaligned="results/{sample}_aligned.bam",
        ref="data/genome.fa"
    output:
        "results/{sample}_merged.bam"
    log:
        "logs/merge_bam_alignment/{sample}.log"
    shell:
        """
        module load biology gatk
        gatk MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ALIGNED_BAM {input.aligned} \
            --UNMAPPED_BAM {input.unaligned} \
            --OUTPUT {wildcards.sample}_merged.bam \
            --REFERENCE_SEQUENCE {input.ref} \
            --PAIRED_RUN true \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --CLIP_ADAPTERS false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --ADD_MATE_CIGAR true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_VERSION "0.7.17-r1188" \
            --PROGRAM_GROUP_COMMAND_LINE "-K 100000000 -p -v 3 -t $SLURM_CPUS_ON_NODE -Y" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true 2> {log}
        """

rule mark_duplicates:
    input:
        "results/{sample}_merged.bam"
    output:
        bam="results/{sample}_marked.bam",
        md5="results/{sample}_marked.bam.md5",
        metrics="results/metrics/{sample}_duplicate_metrics"
    log:
        "logs/mark_duplicates/{sample}.log"
    shell:
        """
        module load biology gatk
        gatk MarkDuplicates \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true 2> {log}
        """

rule sort_and_fix_tags:
    input:
        bam="results/{sample}_marked.bam",
        ref="data/genome.fa"
    output:
        bam="results/{sample}_sorted.bam",
        index="results/{sample}_sorted.bai",
        md5="results/{sample}_sorted.bam.md5"
    log:
        "logs/sort_and_fix_tags/{sample}.log"
    shell:
        """
        set -o pipefail
        module load biology gatk
        gatk SortSam \
            --INPUT {input.bam} \
            --OUTPUT /dev/stdout \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX false \
            --CREATE_MD5_FILE false \
        | \
        gatk SetNmMdAndUqTags \
            --INPUT /dev/stdin \
            --OUTPUT {output.bam} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE {input.ref} 2> {log}
        """

