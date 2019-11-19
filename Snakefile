#SAMPLES = ['A', 'B']
SAMPLES = ['A']

rule all:
    input:
        expand("results/{sample}_merged.bam", sample=SAMPLES)

include: "rules/preprocessing.smk"
