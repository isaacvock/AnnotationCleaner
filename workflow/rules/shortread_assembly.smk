rule SR_stringtie:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        guide=SR_GUIDE_GTF,
    output:
        "results/SR_stringtie/{sample}.gtf",
    log:
        "logs/SR_stringtie/{sample}.log",
    threads: 32
    params:
        extra=ST_STRAND + config["stringtie_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie -G {input.guide} -o {output} -p {threads} {params.extra} {input.bam} 1> {log} 2>&1"



rule SR_stringtie_merge:
    input:
        expand("results/SR_stringtie/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/SR_stringtie_merge/SR_stringtie_merged.gtf",
    log:
        "logs/SR_stringtie_merge/SR_stringtie_merged.log",
    threads: 10
    params:
        extra=config["SR_stringtie_merge_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"

if SR_IS_FINAL:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule SR_filter_annotation:
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

else:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule SR_filter_annotation:
        input:
            "results/clean_assembly/stringtie_merged.gtf",
        output:
            "results/SR_filter_annotation/filtered_SR_annotation.gtf",
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
