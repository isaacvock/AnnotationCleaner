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
        "../scripts/dexseq_prepare_annotation.py 2> {log}"

rule add_exon_reference:
    input:
        "results/raw_flattened_reference/flat_genome.gtf",
    output:
        "results/flattened_reference/flat_genome_exonID.gtf",
    log:
        "logs/add_exon_reference/add_exon.log",
    params:
        shellscript=workflow.source_path("../scripts/exon_ID.sh")
    conda:
        "../envs/add_exon.yaml",
    threads: 1
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {input} {output} 2> {log}
        """


rule quantify_reference_total:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        gtf="results/flattened_reference/flat_genome_exonID.gtf"
    output:
        counts="results/quantify_reference/{sample}_total.csv",
    params:
        strand=config["strandedness"]
    conda:
        "../envs/quantify.yaml"
    log:
        "logs/quantify_reference_total/{sample}.log"
    threads: 1
    shell:
        """
        htseq-count -t aggregate_gene -m intersection-strict -s {params.strand} \
        -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
        -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
        """


rule quantify_reference_exonbin:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        gtf="results/flattened_reference/flat_genome_exonID.gtf"
    output:
        counts="results/quantify_reference/{sample}_exonbin.csv",
    params:
        strand=config["strandedness"]
    conda:
        "../envs/quantify.yaml"
    log:
        "logs/quantify_reference_exonbin/{sample}.log"
    threads: 1
    shell:
        """
        htseq-count -t exonic_part -m union -s {params.strand} \
        -r pos -p bam --add-chromosome-info -i exon_id --nonunique=all \
        -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
        """

rule quantify_reference_exonic:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        gtf="results/flattened_reference/flat_genome_exonID.gtf"
    output:
        counts="results/quantify_reference/{sample}_exonic.csv",
    params:
        strand=config["strandedness"]
    conda:
        "../envs/quantify.yaml"
    log:
        "logs/quantify_reference_exonic/{sample}.log"
    threads: 1
    shell:
        """
        htseq-count -t exonic_part -m intersection-strict -s {params.strand} \
        -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
        -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
        """   

rule clean_reference:
    input:
        ref=config["reference_gtf"],
        cnts_exonic="results/quantify_reference/{sample}_exonic.csv",
        cnts_exonbin="results/quantify_reference/{sample}_exonbin.csv",
        cnts_total="results/quantify_reference/{sample}_total.csv",
    output:
        clean_ref="results/clean_reference/{sample}.gtf"
    params:
        rscript=workflow.source_path("../scripts/clean_reference.R")
    conda:
        "../envs/clean.yaml"
    log:
        "logs/clean_reference/{sample}.log"
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -r {input.ref} -e {input.cnts_exonic} -b {input.cnts_exonbin} -t {input.cnts_total} 1> {log} 2>&1
        """