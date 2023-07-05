rule sort:
    input:
        get_input_bams,
    output:
        "results/sorted/sorted_{sample}.bam"
    log:
        "logs/sorted/{sample}.log"
    threads: 8
    params:
        extra=config["samtools_params"],
    wrapper:
        "v2.1.1/bio/samtools/sort"

rule scallop:
    input:
        "results/sorted/sorted_{sample}.bam",
    output:
        gtf="results/separate_scallops/{sample}.gtf",
    log:
        "logs/separate_scallops/{sample}.log"
    threads: 1
    params:
        extra=config["scallop_params"],
    conda:
        "../envs/scallop.yaml",
    shell:
        "scallop -i {input} -o {output.gtf} {params.extra} 1> {log} 2>&1"
    
rule tacoinput:
    input:
        expand("results/separate_scallops/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/tacoinput/samplefile.txt"
    threads: 1
    run:
        with open(output[0], "w") as file:
            for path in input:
                file.write(path + "\n")
    
rule taco:
    input:
        "results/tacoinput/samplefile.txt"
    output:
        output_dummy="results/ignorethisdirectory/success.txt"
    log:
        "logs/scallop_annotation/scallop.log"
    threads: 24
    params:
        extra=config["taco_params"],
    conda:
        "../envs/scallop.yaml"
    shell:
        """
        taco_run -o ./results/scallop_annotation/ -p {threads} {params.extra} {input} 1> {log} 2>&1
        touch {output.output_dummy}
        """
    