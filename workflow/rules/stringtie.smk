# Use provided reference as guide
if config["use_reference"]:

    rule stringtie:
        input:
            bam="results/sorted/sorted_{sample}.bam",
        output:
            "results/separate_stringties/{sample}.gtf",
        log:
            "logs/separate_stringties/{sample}.log",
        threads: 20
        params:
            extra=ST_STRAND + config["stringtie_params"],
            guide=config["reference_gtf"],
        conda:
            "../envs/stringtie.yaml"
        shell:
            "stringtie -G {params.guide} -o {output} -p {threads} {params.extra} {input.bam} 1> {log} 2>&1"


# Don't use a guide
else:

    rule stringtie:
        input:
            bam="results/sorted/sorted_{sample}.bam",
        output:
            "results/separate_stringties/{sample}.gtf",
        log:
            "logs/separate_stringties/{sample}.log",
        threads: 20
        params:
            extra=ST_STRAND + config["stringtie_params"],
        conda:
            "../envs/stringtie.yaml"
        shell:
            "stringtie -o {output} -p {threads} {params.extra} {input.bam} 1> {log} 2>&1"


rule stringtie_merge:
    input:
        expand("results/separate_stringties/{SID}.gtf", SID=SAMP_NAMES),
    output:
        "results/stringtie_merge/stringtie_merged.gtf",
    log:
        "logs/stringtie_merge/stringtie_merged.log",
    threads: 10
    params:
        extra=config["stringtie_merge_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"


### Remove small transcripts from including SpliceWiz novel targets
rule filter_annotation:
    input:
        "results/clean_assembly/stringtie_merged.gtf",
    output:
        "results/final_annotation/final_annotation.gtf",
    params:
        extra=config["filter_params"],
        rscript=workflow.source_path("../scripts/filter_annotation.R"),
    log:
        "logs/filter_annotation/filter.log",
    conda:
        "../envs/cleaning.yaml"
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -i {input} -o {output} {params.extra} 1> {log} 2>&1
        """
