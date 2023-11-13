# Identify high confidence splice junctions with Portcullis
rule identify_junctions:
    input:
        fasta=config["genome"],
        bams=expand("results/sorted/sorted_{SID}.bam", SID = SAMP_NAMES),
    output:
        "results/identify_junctions/junctions.junctions.bed"
    params:
        extra_prep=config["portcullis_prep_params"],
        extra_junc=config["portcullis_junc_params"],
        strandedness=config["portcullis_strandedness"],
        orientation=config["portcullis_orientation"]
    conda:
        "../envs/portcullis.yaml"
    log:
        "logs/identify_junctions/portcullis.log"
    threads: 8
    shell:
        """
        portcullis prep -t {threads} -o results/prepare_portcullis/ {params.extra_prep} {input.fasta} {input.bams} 1> {log} 2>&1
        portcullis junc -t {threads} --orientation {params.orientation} \portcullis
        --strandedness {params.strandedness} -o results/identify_junctions/junctions {params.extra_junc} \
        results/prepare_portcullis/ 1> {log} 2>&1
        """

# Filter reads associated with low-quality junctions from the bam file
rule filter_bam:
    input:
        junction="results/identify_junctions/junctions.junctions.bed",
        bam="results/sorted/sorted_{sample}.bam",
    output:
        filtered="results/filtered/{sample}_filtered.bam"
    params:
        strand=PORTCULLIS_STRAND,
        clip_mode=config["portcullis_clipping"],
        extra=config["bamfilt_extra"]
    conda:
        "../envs/portcullis.yaml"
    log:
        "logs/filter_bam/{sample}.log"
    threads: 1
    shell:
        "portcullis bamfilt -s {params.strand} -c {params.clip_mode} {params.extra} -o {output.filtered}"