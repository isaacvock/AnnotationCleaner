rule remove_unstranded:
    input:
        gtf="results/stringtie_merge/stringtie_merged.gtf",
    output:
        stranded="results/remove_unstranded/stringtie_merged.gtf",
    log:
        "logs/remove_unstranded/remove_unstranded.log",
    conda:
        "../envs/cleaning.yaml"
    params:
        rscript=workflow.source_path("../scripts/remove_unstranded.R"),
    threads: 1
    shell:
        """
            chmod +x {params.rscript}
            {params.rscript} -i {input.gtf} -o {output.stranded} 1> {log} 2>&1
            """

# DEXSeq flatten merged StringTie assembly
rule flatten_assembly:
    input:
        gtf="results/remove_unstranded/stringtie_merged.gtf",
    output:
        flatgtf="results/raw_flattened_assembly/flat_genome.gtf",
    log:
        "logs/flatten_assembly/flatten.log",
    conda:
        "../envs/dexseq.yaml"
    threads: 1
    script:
        "../scripts/dexseq_prepare_annotation.py"

# Decrease size of exonic bins
rule smaller_bins_assembly:
    input:
        gtf="results/raw_flattened_assembly/flat_genome.gtf",
    output:
        higherres="results/flattened_assembly/flat_genome_binID.gtf",
    log:
        "logs/smaller_bins_assembly/smaller_bins_assembly.log",
    params:
        rscript=workflow.source_path("../scripts/split_bins.R"),
        extra=config["splitting_bins_params"],
    conda:
        "../envs/cleaning.yaml"
    threads: 6
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -i {input.gtf} -o {output.higherres} -t {threads} {params.extra} 1> {log} 2>&1
        """

# Quantify exonic bins
rule quantify_assembly_exonbin:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation="results/flattened_assembly/flat_genome_binID.gtf",
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_assembly/{sample}_exonbin",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 10
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_EB,
    log:
        "logs/quantify_assembly_exonbin/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

# Quantify intronic bins
rule quantify_assembly_intronbin:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation="results/flattened_assembly/flat_genome_binID.gtf",
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_assembly/{sample}_intronbin",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 10
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_IB,
    log:
        "logs/quantify_assembly_intronbin/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

# Clean StringTie assemblies with custom R script
# I am imagining a parameter d that if set, loads sets of csvs as I normally do
rule stringtie_clean_assembly:
    input:
        ref="results/remove_unstranded/stringtie_merged.gtf",
        flatref="results/flattened_assembly/flat_genome_binID.gtf",
        cnts_exonbin=expand(
            "results/quantify_assembly/{SID}_exonbin.featureCounts", SID=SAMP_NAMES
        ),
        cnts_intronbin=expand(
            "results/quantify_assembly/{SID}_intronbin.featureCounts",
            SID=SAMP_NAMES,
        ),
    output:
        clean_ref="results/clean_assembly/stringtie_merged.gtf",
    params:
        rscript=workflow.source_path("../scripts/clean_annotation.R"),
        extra=config["pruning_assembly_params"],
    conda:
        "../envs/cleaning.yaml"
    log:
        "logs/clean_assembly/stringtie_clean_assembly.log",
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -r {input.ref} -f {input.flatref} -d results/quantify_assembly/ -o {output.clean_ref} {params.extra} 1> {log} 2>&1
        """

