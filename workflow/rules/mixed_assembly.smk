rule stringtie_mixed:
    input:
        bams=get_mixed_bams,
        guide=MIX_GUIDE_GTF,
    output:
        "results/stringtie_mixed/{sample}.gtf",
    log:
        "logs/stringtie_mixed/{sample}.log",
    threads: 32
    params:
        extra=ST_STRAND + config["stringtie_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie -G {input.guide} --mix -o {output} -p {threads} {params.extra} {input.bams} 1> {log} 2>&1"


### Merge if using as a guide later
rule stringtie_mixed_merge:
    input:
        expand("results/stringtie_mixed/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/stringtie_mixed_merge/stringtie_mixed_merged.gtf",
    log:
        "logs/stringtie_mixed_merge/stringtie_mixed_merged.log",
    threads: 10
    params:
        extra=config["stringtie_mixed_merge_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"


if MIX_IS_FINAL:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule filter_mixed_annotation:
        input:
            "results/stringtie_mixed_merge/stringtie_mixed_merged.gtf",
        output:
            "results/filter_annotation/filter_mixed_annotation.gtf",
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
    rule filter_mixed_annotation:
        input:
            "results/stringtie_mixed_merge/stringtie_mixed_merged.gtf",
        output:
            "results/filter_mixed_annotation/filter_mixed_annotation.gtf",
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