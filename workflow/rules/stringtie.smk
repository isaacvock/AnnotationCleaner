rule stringtie:
    input:
        "results/sorted/sorted_{sample}.bam",
    output:
        "results/separate_stringties/{sample}.gtf",
    log:
        "logs/separate_stringties/{sample}.log"
    threads: 24
    params:
        extra=config["stringtie_params"],
        gtf=config["reference_gtf"]
    conda:
        "../envs/stringtie.yaml",
    shell:
        "stringtie -o {output} -p {threads} -G {params.gtf} {params.extra} {input} 1> {log} 2>&1"

rule stringtie_merge:
    input:
        expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/stringtie_merge/stringtie_merged.gtf",
    log:
        "logs/stringtie_merge/stringtie_merged.log"
    threads: 24
    params:
        extra=config["stringtie_merge_params"],
        gtf=config["reference_gtf"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge -G {params.gtf} -o {output} {params.extra} {input}"

rule stringtie_tacoinput:
    input:
        expand("results/separate_stringties/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/stringtie_tacoinput/samplefile.txt"
    threads: 1
    run:
        with open(output[0], "w") as file:
            for path in input:
                file.write(path + "\n")

rule stringtie_taco:
    input:
        "results/stringtie_tacoinput/samplefile.txt"
    output:
        output_dummy="results/ignorethisdirectory_stringtie/success.txt"
    log:
        "logs/stringtie_taco/stringtie_taco.log"
    threads: 24
    params:
        extra_taco=config["taco_params"],
        extra_refcomp=config["refcomp_params"],
        gtf=config["reference_gtf"]
    conda:
        "../envs/scallop.yaml"
    shell:
        """
        taco_run -o ./results/stringtie_taco/ -p {threads} {params.extra_taco} {input} 1> {log} 2>&1
        taco_refcomp -o ./results/stringtie_taco_refcomp/ -r {params.gtf} -t ./results/stringtie_taco/assembly.gtf {params.extra_refcomp} 1> {log} 2>&1
        touch {output.output_dummy}
        """
    