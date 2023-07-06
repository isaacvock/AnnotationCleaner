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
    
rule scallop_tacoinput:
    input:
        expand("results/separate_scallops/{SID}.gtf", SID = SAMP_NAMES),
    output:
        "results/scallop_tacoinput/samplefile.txt"
    threads: 1
    run:
        with open(output[0], "w") as file:
            for path in input:
                file.write(path + "\n")
    
rule scallop_taco:
    input:
        "results/scallop_tacoinput/samplefile.txt"
    output:
        output_dummy="results/ignorethisdirectory/success.txt"
    log:
        "logs/scallop_taco/scallop.log"
    threads: 24
    params:
        extra_taco=config["taco_params"],
        extra_refcomp=config["refcomp_params"],
        gtf=config["reference_gtf"]
    conda:
        "../envs/scallop.yaml"
    shell:
        """
        taco_run -o ./results/scallop_taco/ -p {threads} {params.extra_taco} {input} 1> {log} 2>&1
        taco_refcomp -o ./results/scallop_taco_refcomp/ -r {params.gtf} -t ./results/scallop_taco/assembly.gtf {params.extra_refcomp} 1> {log} 2>&1
        touch {output.output_dummy}
        """
    