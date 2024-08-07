import glob

# Sample names, used in multiple places to help Snakemake infer wildcards
# and to expand list of bam files
SAMP_NAMES = list(config['samples'].keys())

# Retrieve input bam files for first steps
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]

# List of paths to Stringtie outputs for Mikado
if config["stringtie"]["use_stringtie"]:
    STRINGTIE_PATHS = expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES)    

# Target rule (so final output to be looked for)
def get_target_input():

    target = []

    if config["clean_only"]:

        target.append("results/clean_reference/cleaned_reference.gtf")

    else:

        if config["stringtie"]["use_stringtie"]:

            target.append(expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES))

        if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_taco"]:

            target.append("results/ignorethisdirectory_stringtie/success.txt")

        if config["stringtie"]["use_stringtie"] and config["stringtie"]["use_merge"]:
        
            target.append("results/stringtie_merge/stringtie_merged.gtf")

        # Filtered annotation
        target.append("results/final_annotation/final_annotation.gtf")

    return target


# # IDs for all split up sub-fasta files
# num_digits = len(str(config["num_sub"]))

# SPLIT_IDS = [str(i).zfill(num_digits) for i in range(0, config["num_sub"])]


# # Number of threads that can be used for serialization step
# if config["num_sub"] > 1:
#     SERIALISE_THREADS = config["num_sub"]

#     if SERIALISE_THREADS > workflow.cores:
#         SERIALISE_THREADS = workflow.cores

# else:
#     SERIALISE_THREADS = 1


### FeatureCounts helpers

# Strandedness inference
if config["strandedness"] == "yes":

    STRANDEDNESS = 1

elif config["strandedness"] == "reverse":

    STRANDEDNESS = 2

else:

    STRANDEDNESS = 0


# Extra parameters

FC_EXTRA = config["feature_counts_params"]

if config["PE"]:

    FC_EXTRA = FC_EXTRA + " -O -f -p"

else:

    FC_EXTRA = FC_EXTRA + " -O -f"



FC_EXTRA_IB = FC_EXTRA + " -g intron_id -t intronic_part"
FC_EXTRA_EB = FC_EXTRA + " -g exon_id -t exonic_part"





### Sort bam files
rule sort:
    input:
        get_input_bams,
    output:
        "results/sorted/sorted_{sample}.bam"
    log:
        "logs/sort/{sample}.log"
    threads: 8
    params:
        extra=config["samtools_params"],
    wrapper:
        "v2.1.1/bio/samtools/sort"