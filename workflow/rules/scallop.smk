rule sort:
    input:
        get_input_bams,
    output:
        "results/sorted/sorted_{sample}.bam"
    log:
        "logs/sorted/{sample}.log"
    threads: 8
    params:
        extra=config["samtools_params],
    wrapper:
        "v2.1.1/bio/samtools/sort"

rule scallop:
    input:
        "results/sorted/sorted_{sample}.bam",
    output:
        "results/separate_scallops/{sample}.gtf"
    log:
        "logs/separate_scallops/{sample}.log"
    threads: 1
    params:
        extra=config["scallop_params"],
    conda:
        "../envs/scallop.yaml",
    shell:
        "scallop -i {input} -o {output} {params.extra} 1> {log} 2>&1"
    
rule taco:
    input:
        expand("results/separate_scallops/{SID}.gtf", SID = SAMP_NAMES)
    output:
        directory("results/scallop_annotation/")
    log:
        "logs/scallop_annotation/scallop.log"
    threads: 24
    params:
        extra=config["taco_params"],
    conda:
        "../envs/scallop.yaml"
    shell:
        "taco_run {input} -o {output} -p {threads} {params.extra} 1> {log} 2>&1"
    