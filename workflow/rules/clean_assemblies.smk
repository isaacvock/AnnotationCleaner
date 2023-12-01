rule remove_unstranded:
    input:
        gtf="results/separate_stringties/{sample}.gtf"
    output:
        stranded="results/remove_unstranded/{sample}.gtf"
    log:
        "logs/remove_unstranded/{sample}.log",
    conda:
        "../envs/cleaning.yaml"
    params:
        rscript=workflow.source_path("../scripts/remove_unstranded.R")
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -i {input.gtf} -o {output.stranded} 1> {log} 2>&1
        """


if config["stringtie"]["clean_then_merge"]:

    # DEXSeq flatten individual StringTie assemblies
    rule flatten_assembly:
        input:
            gtf="results/remove_unstranded/{sample}.gtf"
        output:
            flatgtf="results/raw_flattened_assembly/{sample}_flat_genome.gtf",
        log:
            "logs/flatten_assembly/{sample}.log",
        conda: 
            "../envs/dexseq.yaml"
        threads: 1
        script:
            "../scripts/dexseq_prepare_annotation.py"

    # Decrease size of exonic bins
    rule smaller_bins_assembly:
        input:
            gtf="results/raw_flattened_assembly/{sample}_flat_genome.gtf",
        output:
            higherres="results/flattened_assembly/{sample}_flat_genome_binID.gtf",
        log:
            "logs/smaller_bins_assembly/{sample}.log"
        params:
            rscript=workflow.source_path("../scripts/split_bins.R"),
            extra=config["splitting_bins_params"]
        conda:
            "../envs/cleaning.yaml"
        threads: 6
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -i {input.gtf} -o {output.higherres} -t {threads} {params.extra} 1> {log} 2>&1
            """

    # Total read count quantification for each gene in StringTie assemblies
    rule quantify_assembly_total:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/{sample}_flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_total.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_total/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t aggregate_gene -m intersection-strict -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Exon bin quantification for each gene in StringTie assemblies
    rule quantify_assembly_exonbin:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/{sample}_flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_exonbin.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_exonbin/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t exonic_part -m union -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i exon_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Intron bin quantification for each gene in StringTie assemblies
    rule quantify_assembly_intronbin:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/{sample}_flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_intronbin.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_intronbin/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t intronic_part -m union -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i intron_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Exonic quantification for each gene in StringTie assemblies
    rule quantify_assembly_exonic:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/{sample}_flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_exonic.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_exonic/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t exonic_part -m intersection-strict -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """ 

    # Clean StringTie assemblies with custom R script
    rule stringtie_clean_assembly:
        input:
            ref="results/remove_unstranded/{sample}.gtf",
            flatref="results/flattened_assembly/{sample}_flat_genome_binID.gtf",
            cnts_exonic="results/quantify_assembly/{sample}_exonic.csv",
            cnts_exonbin="results/quantify_assembly/{sample}_exonbin.csv",
            cnts_total="results/quantify_assembly/{sample}_total.csv",
            cnts_intronbin="results/quantify_assembly/{sample}_intronbin.csv",
        output:
            clean_ref="results/clean_assembly/{sample}.gtf"
        params:
            rscript=workflow.source_path("../scripts/clean_annotation.R"),
            extra=config["pruning_assembly_params"]
        conda:
            "../envs/cleaning.yaml"
        log:
            "logs/clean_assembly/{sample}.log"
        threads: 1
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -r {input.ref} -f {input.flatref} -e {input.cnts_exonic} -b {input.cnts_exonbin} -t {input.cnts_total} \
            -d ./results/quantify_assembly/ -u {output.cnts_intronbin} -o {output.clean_ref} {params.extra} 1> {log} 2>&1
            """

else:

    # DEXSeq flatten merged StringTie assembly
    rule flatten_assembly:
        input:
            gtf="results/stringtie_merge/stringtie_merged.gtf"
        output:
            flatgtf="results/raw_flattened_assembly/flat_genome.gtf",
        log:
            "logs/flatten_assembly/flatten.log",
        conda: 
            "../envs/dexseq.yaml"
        threads: 1
        script:
            "../scripts/dexseq_prepare_annotation.py"

    # Decrease size of exonic bins
    rule smaller_bins_assembly:
        input:
            gtf="results/raw_flattened_assembly/flat_genome.gtf",
        output:
            higherres="results/flattened_assembly/flat_genome_binID.gtf",
        log:
            "logs/smaller_bins_assembly/{sample}.log"
        params:
            rscript=workflow.source_path("../scripts/split_bins.R"),
            extra=config["splitting_bins_params"]
        conda:
            "../envs/cleaning.yaml"
        threads: 6
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -i {input.gtf} -o {output.higherres} -t {threads} {params.extra} 1> {log} 2>&1
            """

    # Total read count quantification for each gene in StringTie assemblies
    rule quantify_assembly_total:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_total.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_total/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t aggregate_gene -m intersection-strict -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Exon bin quantification for each gene in StringTie assemblies
    rule quantify_assembly_exonbin:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_exonbin.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_exonbin/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t exonic_part -m union -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i exon_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Intron bin quantification for each gene in StringTie assemblies
    rule quantify_assembly_intronbin:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_intronbin.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_exonbin/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t intronic_part -m union -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i intron_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """

    # Exonic quantification for each gene in StringTie assemblies
    rule quantify_assembly_exonic:
        input:
            bam="results/sorted/sorted_{sample}.bam",
            gtf="results/flattened_assembly/flat_genome_binID.gtf"
        output:
            counts="results/quantify_assembly/{sample}_exonic.csv",
        params:
            strand=config["strandedness"]
        conda:
            "../envs/quantify.yaml"
        log:
            "logs/quantify_assembly_exonic/{sample}.log"
        threads: 1
        shell:
            """
            htseq-count -t exonic_part -m intersection-strict -s {params.strand} \
            -r pos -p bam --add-chromosome-info -i gene_id --nonunique=all \
            -c {output.counts} {input.bam} {input.gtf} 1> {log} 2>&1
            """ 

    # Clean StringTie assemblies with custom R script
        # I am imagining a parameter d that if set, loads sets of csvs as I normally do
    rule stringtie_clean_assembly:
        input:
            ref="results/remove_unstranded/{sample}.gtf",
            flatref="results/flattened_assembly/flat_genome_binID.gtf",
            cnts_exonic=expand("results/quantify_assembly/{SID}_exonic.csv", SID = SAMP_NAMES),
            cnts_exonbin=expand("results/quantify_assembly/{SID}_exonbin.csv", SID = SAMP_NAMES),
            cnts_total=expand("results/quantify_assembly/{SID}_total.csv", SID = SAMP_NAMES),
            cnts_intronbin=expand("results/quantify_assembly/{SID}_intronbin.csv", SID = SAMP_NAMES)
        output:
            clean_ref="results/clean_assembly/stringtie_merged.gtf"
        params:
            rscript=workflow.source_path("../scripts/clean_annotation.R"),
            extra=config["pruning_assembly_params"]
        conda:
            "../envs/cleaning.yaml"
        log:
            "logs/clean_assembly/{sample}.log"
        threads: 1
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -r {input.ref} -f {input.flatref} -d results/quantify_assembly/ -o {output.clean_ref} {params.extra} 1> {log} 2>&1
            """