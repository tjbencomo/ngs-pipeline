SAMPLES = ['A', 'B']

rule all:
    input:
        expand("results/{sample}_aligned.bam", sample=SAMPLES)

include: "rules/preprocessing.smk"
