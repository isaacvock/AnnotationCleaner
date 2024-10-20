rule mixed_stringtie:
    input:
        bams=get_mixed_bams,
        guide=MIX_GUIDE_GTF,
    output:
        "results/mixed_stringtie/{sample}.gtf",
    log:
        "logs/mixed_stringtie/{sample}.log",
    threads: 32
    params:
        extra=ST_STRAND + config["mixed_stringtie_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie -G {input.guide} --mix -o {output} -p {threads} {params.extra} {input.bams} 1> {log} 2>&1"


### Merge if using as a guide later
rule mixed_stringtie_merge:
    input:
        expand("results/mixed_stringtie/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/mixed_stringtie_merge/mixed_stringtie_merged.gtf",
    log:
        "logs/mixed_stringtie_merge/mixed_stringtie_merged.log",
    threads: 10
    params:
        extra=config["mixed_stringtie_merge_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"


if MIX_IS_FINAL:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule mixed_filter_annotation:
        input:
            FINAL_INPUT
        output:
            "results/final_annotation/final_annotation.gtf",
        params:
            extra=config["mixed_filter_params"],
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
    rule mixed_filter_annotation:
        input:
            "results/mixed_stringtie_merge/mixed_stringtie_merged.gtf",
        output:
            "results/mixed_filter_annotation/mixed_filter_annotation.gtf",
        params:
            extra=config["mixed_filter_params"],
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