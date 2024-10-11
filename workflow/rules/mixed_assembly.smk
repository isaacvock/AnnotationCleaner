if config["use_reference"] or LR_GUIDED:

    rule stringtie_mixed:
        input:
            bams=get_mixed_bams,
            guide=GUIDE_GTF,
        output:
            "results/stringtie_mixed/{sample}.gtf",
        log:
            "logs/stringtie_mixed/{sample}.log",
        threads: 20
        params:
            extra=ST_STRAND + config["stringtie_params"],
        conda:
            "../envs/stringtie.yaml"
        shell:
            "stringtie -G {input.guide} --mix -o {output} -p {threads} {params.extra} {input.bams} 1> {log} 2>&1"

else:

    rule stringtie_mixed:
        input:
            get_mixed_bams
        output:
            "results/stringtie_mixed/{sample}.gtf",
        log:
            "logs/stringtie_mixed/{sample}.log",
        threads: 20
        params:
            extra=ST_STRAND + config["stringtie_params"],
        conda:
            "../envs/stringtie.yaml"
        shell:
            "stringtie --mix -o {output} -p {threads} {params.extra} {input} 1> {log} 2>&1"
