import glob

SAMP_NAMES = list(config['samples'].keys())

def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]

def get_mikado_input():

    gtfs = []

    if config["scallop"]["use_scallop"]:
        gtfs.append("results/ignorethisdirectory_scallop/success.txt")
    
    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_taco"]:
        gtfs.append("results/ignorethisdirectory_stringtie/success.txt")

    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_merge"]:
        gtfs.append("results/stringtie_merge/stringtie_merged.gtf")


    return gtfs


def get_target_input():

    target = []

    if config["use_mikado"]:

        target.append("results/ignorethisdirectory_mikado/success.txt")

    if config["scallop"]["use_scallop"]:

        target.append("results/ignorethisdirectory_scallop/success.txt")

    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_taco"]:

        target.append("results/ignorethisdirectory_stringtie/success.txt")

    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_merge"]:
    
        target.append("results/stringtie_merge/stringtie_merged.gtf")

    return target



### Sort bam files
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