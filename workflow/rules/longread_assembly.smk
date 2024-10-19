rule longread_stringtie:
    input:
        bam="results/sorted/sorted_{sample}.bam",
        guide=LR_GUIDE_GTF
    output:
        "results/longread_stringtie/{sample}.gtf",
    log:
        "logs/longread_stringtie/{sample}.log",
    threads: 32
    params:
        extra=config["stringtie_longread_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie -G {input.guide} -L -o {output} -p {threads} {params.extra} {input.bam} 1> {log} 2>&1"


rule longread_stringtie_merge:
    input:
        expand("results/longread_stringtie/{SID}.gtf", SID=LONGREAD_NAMES),
    output:
        "results/longread_stringtie_merge/stringtie_merged.gtf",
    log:
        "logs/longread_stringtie_merge/stringtie_merged.log",
    threads: 10
    params:
        extra=config["stringtie_merge_longread_params"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -p {threads} -o {output} {params.extra} {input}"


if LR_IS_FINAL:

    ### Remove small transcripts from including SpliceWiz novel targets
    rule filter_longread_annotation:
        input:
            "results/longread_stringtie_merge/stringtie_merged.gtf",
        output:
            "results/final_annotation/final_annotation.gtf",
        params:
            extra=config["filter_params"],
            rscript=workflow.source_path("../scripts/filter_annotation.R"),
        log:
            "logs/filter_longread_annotation/filter.log",
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
    rule filter_longread_annotation:
        input:
            "results/longread_stringtie_merge/stringtie_merged.gtf",
        output:
            "results/longread_stringtie_merge/filtered_longread_annotation.gtf",
        params:
            extra=config["filter_params"],
            rscript=workflow.source_path("../scripts/filter_annotation.R"),
        log:
            "logs/filter_longread_annotation/filter.log",
        conda:
            "../envs/cleaning.yaml"
        threads: 1
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -i {input} -o {output} {params.extra} 1> {log} 2>&1
            """