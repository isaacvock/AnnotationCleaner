import glob

# Sample names, used in multiple places to help Snakemake infer wildcards
# and to expand list of bam files
SAMP_NAMES = list(config['samples'].keys())

# Retrieve input bam files for first steps
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]

# List of paths to Scallop outputs for Mikado
if config["scallop"]["use_scallop"]:
    SCALLOP_PATHS = expand("results/separate_scallops/{SID}.gtf", SID = SAMP_NAMES)

# List of paths to Stringtie outputs for Mikado
if config["stringtie"]["use_stringtie"]:
    STRINGTIE_PATHS = expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES)    

# Retrieve lists of files that must get created for mikado to run
def get_mikado_input():

    gtfs = []

    if config["scallop"]["use_scallop"]:
        gtfs.append(expand("results/separate_scallops/{SID}.gtf", SID = SAMP_NAMES))
    
    if config["stringtie"]["use_stringtie"]:
        gtfs.append(expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES))

    return gtfs

# Target rule (so final output to be looked for)
def get_target_input():

    target = []

    if config["use_mikado"]:

        target.append("results/mikado_pick/mikado.loci.gff3")

    if config["scallop"]["use_scallop"] and config["scallop"]["use_taco"]:

        target.append("results/ignorethisdirectory_scallop/success.txt")

    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_taco"]:

        target.append("results/ignorethisdirectory_stringtie/success.txt")

    if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_merge"]:
    
        target.append("results/stringtie_merge/stringtie_merged.gtf")

    return target


# IDs for all split up sub-fasta files
SPLIT_IDS =  list(range(1, config["num_sub"] + 1))


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