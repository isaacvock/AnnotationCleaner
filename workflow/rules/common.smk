import glob

SAMP_NAMES = list(config['samples'].keys())

def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]

rule sort:
    input:
        get_input_bams,
    output:
        "results/sorted/sorted_{sample}.bam"
    log:
        "logs/sorted/{sample}.log"
    threads: 8
    params:
        extra=config["samtools_params"],
    wrapper:
        "v2.1.1/bio/samtools/sort"