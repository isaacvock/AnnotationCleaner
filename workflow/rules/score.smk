rule quantify_finalGTF_exonic:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation=score_input,
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_intronic_coverage/{sample}_exonic",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_EXON,
    log:
        "logs/quantify_finalGTF_exonic/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


rule quantify_finalGTF_gene:
    input:
        # list of sam or bam files
        samples="results/sorted/sorted_{sample}.bam",
        annotation=score_input,
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/quantify_intronic_coverage/{sample}_gene",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 20
    params:
        strand=STRANDEDNESS,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=FC_EXTRA_GENE,
    log:
        "logs/quantify_finalGTF_gene/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

