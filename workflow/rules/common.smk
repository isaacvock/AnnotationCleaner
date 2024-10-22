import glob

# Sample names, used in multiple places to help Snakemake infer wildcards
# and to expand list of bam files
SAMP_NAMES = list(config["samples"].keys())

# Were long reads provided?
LONGREADS_PROVIDED = len(config["long_reads"]) > 0
LONGREAD_NAMES = []

if LONGREADS_PROVIDED:

    LONGREAD_NAMES = list(set(SAMP_NAMES) & set(config["long_reads"]))
    SAMP_NAMES = list(set(SAMP_NAMES) - set(config["long_reads"]))


### INFER PIPELINE PATHWAY FROM CONFIG VALUES

# Possible paths for various inputs/outputs
DUMMY_PATH = "results/dummypath"
REFERENCE = config["reference_gtf"]
CLEANED_REFERENCE = "results/clean_reference/cleaned_reference.gtf"
CLEANED_ASSEMBLY = "results/clean_assembly/cleaned_assembly.gtf"
FILTERED_LR_GTF = "results/longread_filter_annotation/longread_filtered_annotation.gtf"
FILTERED_SR_GTF = "results/shortread_filter_annotation/shortread_filtered_annotation.gtf"
FILTERED_MIX_GTF = "results/mixed_filter_annotation/mixed_filter_annotation.gtf"
MERGED_MIX_GTF = "results/mixed_stringtie_merge/mixed_stringtie_merged.gtf"
MERGED_SR_GTF = "results/shortread_stringtie_merge/shortread_stringtie_merged.gtf"

# Pathway inference
LRSR_STRAT = config["options"]["LRSR_strat"]

if LRSR_STRAT == "mix_only":

    config["LRonly_first"] = False
    config["use_mix"] = True
    MIX_IS_FINAL = True
    SR_IS_FINAL = False
    LR_IS_FINAL = False
    LR_GUIDE_GTF = DUMMY_PATH
    SR_GUIDE_GTF = DUMMY_PATH

    if config["options"]["trim_reference"]:

        MIX_GUIDE_GTF = CLEANED_REFERENCE

    else:

        MIX_GUIDE_GTF = REFERENCE

    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_MIX_GTF

    else:

        FINAL_INPUT = MERGED_MIX_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH


elif LRSR_STRAT == "LR_then_SR":

    config["LRonly_first"] = False
    config["use_mix"] = False
    MIX_IS_FINAL = False
    SR_IS_FINAL = True
    LR_IS_FINAL = False

    if config["options"]["trim_reference"]:

        LR_GUIDE_GTF = CLEANED_REFERENCE


    else:

        LR_GUIDE_GTF = REFERENCE

    SR_GUIDE_GTF = FILTERED_LR_GTF
    MIX_GUIDE_GTF = DUMMY_PATH

    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_SR_GTF

    else:

        FINAL_INPUT = MERGED_SR_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH

elif LRSR_STRAT == "mix_then_SR":

    config["LRonly_first"] = False
    config["use_mix"] = True
    MIX_IS_FINAL = False
    SR_IS_FINAL = True
    LR_IS_FINAL = False
    LR_GUIDE_GTF = DUMMY_PATH
    SR_GUIDE_GTF = FILTERED_MIX_GTF

    if config["options"]["trim_reference"]:

        MIX_GUIDE_GTF = CLEANED_REFERENCE

    else:

        MIX_GUIDE_GTF = REFERENCE


    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_SR_GTF

    else:

        FINAL_INPUT = MERGED_SR_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH


elif LRSR_STRAT == "LR_then_mix":

    config["LRonly_first"] = True
    config["use_mix"] = True
    MIX_IS_FINAL = True
    SR_IS_FINAL = False
    LR_IS_FINAL = False
    
    if config["options"]["trim_reference"]:

        LR_GUIDE_GTF = CLEANED_REFERENCE

    else:

        LR_GUIDE_GTF = REFERENCE

    SR_GUIDE_GTF = DUMMY_PATH
    MIX_GUIDE_GTF = FILTERED_LR_GTF


    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_MIX_GTF

    else:

        FINAL_INPUT = MERGED_MIX_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH

elif LRSR_STRAT == "SR_then_mix":

    config["LRonly_first"] = False
    config["use_mix"] = True
    MIX_IS_FINAL = True
    SR_IS_FINAL = False
    LR_IS_FINAL = False
    
    if config["options"]["trim_reference"]:

        SR_GUIDE_GTF = CLEANED_REFERENCE

    else:

        SR_GUIDE_GTF = REFERENCE

    LR_GUIDE_GTF = DUMMY_PATH
    MIX_GUIDE_GTF = FILTERED_SR_GTF


    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_MIX_GTF

    else:

        FINAL_INPUT = MERGED_MIX_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH

# If long read data isn't provided:
if not LONGREADS_PROVIDED:

    config["LRonly_first"] = False
    config["use_mix"] = False
    MIX_IS_FINAL = False
    SR_IS_FINAL = True
    LR_IS_FINAL = False

    if config["options"]["trim_reference"]:

        SR_GUIDE_GTF = CLEANED_REFERENCE

    else:

        SR_GUIDE_GTF = REFERENCE


    if config["options"]["trim_assembly"]:

        FINAL_INPUT = CLEANED_ASSEMBLY
        ASSEMBLY_CLEANING_INPUT = MERGED_SR_GTF

    else:

        FINAL_INPUT = MERGED_SR_GTF
        ASSEMBLY_CLEANING_INPUT = DUMMY_PATH


### END OF PIPELINE PATHWAY INFERENCE


# Retrieve input bam files for first steps
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# Get bam file combo for longread + short read mixed assembly
def get_mixed_bams(wildcards):
    LR_str = str(config["LRSR_pairs"][wildcards.sample])

    SR_bam = f"results/sorted/sorted_{wildcards.sample}.bam"
    LR_bam = f"results/sorted/sorted_{LR_str}.bam"

    return [SR_bam, LR_bam]



# Target rule (so final output to be looked for)
def get_target_input():
    target = []

    if config["clean_only"]:
        target.append("results/clean_reference/cleaned_reference.gtf")

    else:
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
    shell:
        """
        samtools sort -@ {threads} {params.extra} -o {output} {input} 1> {log} 2>&1
        """


### For intron fraction calculation

if config["clean_only"]:
    score_input = "results/clean_reference/cleaned_reference.gtf"

else:
    score_input = "results/final_annotation/final_annotation.gtf"
