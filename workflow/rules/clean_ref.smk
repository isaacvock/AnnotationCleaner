rule flatten_reference:
    input:
        gtf=config["reference_gtf"]
    output:
        flatgtf="results/raw_flattened_reference/flat_genome.gtf",
    log:
        "logs/flatten_reference/flatten.log",
    conda: 
        "../envs/dexseq.yaml"
    threads: 1
    script:
        "../scripts/dexseq_prepare_annotation.py"

rule smaller_bins_reference:
    input:
        gtf="results/raw_flattened_reference/flat_genome.gtf",
    output:
        higherres=config["flat_ref"],
    log:
        "logs/smaller_bins_reference/smaller_bins.log"
    params:
        rscript=workflow.source_path("../scripts/split_bins.R"),
        extra=config["splitting_bins_params"]
    conda:
        "../envs/cleaning.yaml"
    threads: 6
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -i {input.gtf} -o {output.higherres} -t {threads} {params.extra} 1> {log} 2>&1
        """


# Quantify exonic bins
rule quantify_reference_exonbin:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation=config["flat_ref"],
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_reference/{sample}_exonbin",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 10
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_EB,
    log:
        "logs/quantify_reference_exonbin/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

# Quantify intronic bins
rule quantify_reference_intronbin:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation=config["flat_ref"],
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_reference/{sample}_intronbin",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 10
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_IB,
    log:
        "logs/quantify_reference_intronbin/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


if config["clean_only"]:

    rule clean_reference:
        input:
            ref=config["reference_gtf"],
            flatref=config["flat_ref"],
            cnts_exonbin=expand("results/quantify_reference/{SID}_exonbin.featureCounts", SID = SAMP_NAMES),
            cnts_intronbin=expand("results/quantify_reference/{SID}_intronbin.featureCounts", SID = SAMP_NAMES),
        output:
            clean_ref="results/clean_reference/cleaned_reference.gtf"
        params:
            rscript=workflow.source_path("../scripts/clean_annotation.R"),
            extra=config["pruning_reference_params"]
        conda:
            "../envs/cleaning.yaml"
        log:
            "logs/clean_reference/clean_reference.log"
        threads: 1
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -r {input.ref} -f {input.flatref} -b {input.cnts_exonbin} \
            -u {input.cnts_intronbin} -o {output.clean_ref} -d ./results/quantify_reference/ {params.extra} 1> {log} 2>&1
            """

else:

    rule clean_reference:
        input:
            ref=config["reference_gtf"],
            flatref=config["flat_ref"],
            cnts_exonbin="results/quantify_reference/{sample}_exonbin.featureCounts",
            cnts_intronbin="results/quantify_reference/{sample}_intronbin.featureCounts"
        output:
            clean_ref="results/clean_reference/{sample}.gtf"
        params:
            rscript=workflow.source_path("../scripts/clean_annotation.R"),
            extra=config["pruning_reference_params"]
        conda:
            "../envs/cleaning.yaml"
        log:
            "logs/clean_reference/{sample}.log"
        threads: 1
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -r {input.ref} -f {input.flatref} -b {input.cnts_exonbin} \
            -u {input.cnts_intronbin} -o {output.clean_ref} -d ./results/quantify_reference/ {params.extra} 1> {log} 2>&1
            """