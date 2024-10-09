import glob

# Sample names, used in multiple places to help Snakemake infer wildcards
# and to expand list of bam files
SAMP_NAMES = list(config["samples"].keys())


# Retrieve input bam files for first steps
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# List of paths to Stringtie outputs for Mikado
if not config["clean_only"]:
    STRINGTIE_PATHS = expand("results/separate_stringties/{SID}.gtf", SID=SAMP_NAMES)


# Target rule (so final output to be looked for)
def get_target_input():
    target = []

    if config["clean_only"]:
        target.append("results/clean_reference/cleaned_reference.gtf")

    else:
        # Each StringTie assembly
        target.append(expand("results/separate_stringties/{SID}.gtf", SID=SAMP_NAMES))

        # Merged StringTie assembly
        target.append("results/stringtie_merge/stringtie_merged.gtf")

        # Filtered annotation
        target.append("results/final_annotation/final_annotation.gtf")

    # For scoring final annotations
    target.append(
        expand(
            "results/quantify_intronic_coverage/{SID}_gene.featureCounts",
            SID=SAMP_NAMES,
        )
    )
    target.append(
        expand(
            "results/quantify_intronic_coverage/{SID}_exonic.featureCounts",
            SID=SAMP_NAMES,
        )
    )

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


### StringTie helpers
if config["strandedness"] == "yes":
    ST_STRAND = "--fr "

elif config["strandedness"] == "reverse":
    ST_STRAND = "--rf "

else:
    ST_STRAND = ""


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
    FC_EXTRA = FC_EXTRA + " -O -p"

else:
    FC_EXTRA = FC_EXTRA + " -O"


FC_EXTRA_IB = FC_EXTRA + " -f -g intron_id -t intronic_part --extraAttributes gene_id"
FC_EXTRA_EB = FC_EXTRA + " -f -g exon_id -t exonic_part --extraAttributes gene_id"
FC_EXTRA_GENE = FC_EXTRA + " -g gene_id -t transcript"

NONOVERLAP = config["feature_counts_exon_nonoverlap"]
FC_EXTRA_EXON = FC_EXTRA + " -g gene_id --nonOverlap " + NONOVERLAP


### Sort bam files
rule sort:
    input:
        get_input_bams,
    output:
        "results/sorted/sorted_{sample}.bam",
    log:
        "logs/sort/{sample}.log",
    threads: 8
    params:
        extra=config["samtools_params"],
    conda:
        "../envs/sort.yaml"
    script:
        "../scripts/samtools_sort.py"


### For intron fraction calculation

if config["clean_only"]:
    score_input = "results/clean_reference/cleaned_reference.gtf"

else:
    score_input = "results/final_annotation/final_annotation.gtf"
