rule shortread_stringtie:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        guide=SR_GUIDE_GTF,
    output:
        "results/shortread_stringtie/{sample}.gtf",
    log:
        "logs/shortread_stringtie/{sample}.log",
    threads: 32
    params:
        extra=ST_STRAND + config["shortread_stringtie_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie -G {input.guide} -o {output} -p {threads} {params.extra} {input.bam} 1> {log} 2>&1"



rule shortread_stringtie_merge:
    input:
        expand("results/SR_stringtie/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/shortread_stringtie_merge/shortread_stringtie_merged.gtf",
    log:
        "logs/shortread_stringtie_merge/shortread_stringtie_merge.log",
    threads: 10
    params:
        extra=config["shortread_stringtie_merge_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"

if SR_IS_FINAL:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule shortread_filter_annotation:
        input:
            FINAL_INPUT
        output:
            "results/final_annotation/final_annotation.gtf",
        params:
            extra=config["shortread_filter_params"],
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
    rule shortread_filter_annotation:
        input:
            "results/shortread_stringtie_merge/shortead_stringtie_merged.gtf",
        output:
            "results/shortread_filter_annotation/shortread_filtered_annotation.gtf",
        params:
            extra=config["shortread_filter_params"],
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
