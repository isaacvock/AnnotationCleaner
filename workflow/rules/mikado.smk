# Identify high confidence splice junctions with Portcullis
rule identify_junctions:
    input:
        fasta=config["genome"],
        bams=expand("results/sorted/sorted_{SID}.bam", SID = SAMP_NAMES)
    output:
        "results/identify_junctions/junctions.bed"
    params:
        extra_prep=config["portcullis_prep_params"],
        extra_junc=config["portcullis_junc_params"],
        strandedness=config["portcullis_strandedness"],
        orientation=config["portcullis_orientation"]
    conda:
        "../envs/mikado.yaml"
    threads: 8
    shell:
        """
        portcullis prep -t {threads} -o prepare_portcullis {params.extra_prep} {input.fasta} {input.bams}
        portcullis junc -t {threads} --orientation {params.orientation} \
        --strandedness {params.strandedness} -o results/identify_junctions/junctions {params.extra_junc} \
        prepare_portcullis/
        """

# Create configuration file
rule mikado_configure:
    input:
        mlist=config["mikado_list"],
        reference=config["genome"],
        junctions="results/identify_junctions/junctions.bed"
    output:
        "results/mikado_configure/configuration.yaml"
    params:
        extra=config["mikado_configure_params"]
    threads: 4
    conda:
        "../envs/mikado.yaml"
    shell:
        "mikado configure --list {input.mlist} --reference {input.reference} --junctions {input.junctions} -od results/mikado_configure/ {params.extra} -t {threads}"  

# Create sorted, non-redundant GTF with all input assemblies
rule mikado_prepare:
    input:
        "results/mikado_configure/configuration.yaml",
    output:
        gtf="results/mikado_prepare/mikado_prepared.gtf",
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    params:
        extra=config["mikado_prepare_params"]
    conda:
        "../envs/mikado.yaml"
    shell:
        "mikado prepare --json-conf {input} -od results/mikado_prepare/ {params.extra}"

rule identify_orfs:
    input:
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    output:
        orfs="results/identify_orfs/mikado.orfs.gff3"
    params:
        orf_length=config["orf_min_length"],
        extra=config["transdecoder_params"],
    conda:
        "../envs/mikado.yaml"
    shell:
        "TransDecoder.LongOrfs -t {input.fasta} -m {params.orf_length} --output_dir results/identify_orfs/ {params.extra}"


# Run BLAST to get homology data that will help mikado
rule mikado_blast:
    input:
        proteins=config["protein_fasta"],
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    output:
        prepare_log="results/mikado_blast/blast_prepare.log",
        mikado_blast="results/mikado_blast/mikado_prepared.blast.tsv",
    conda:
        "../envs/mikado.yaml"
    shell:
        """
        makeblastdb -in {input.proteins} -dbtype prot -parse_seqids > {output.prepare_log}
        blastx -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qent sstart send evalue bitscore ppos btop" \
        -num_threads {threads} -query mikado_prepared.fasta -db {input.fasta} -out {output.mikado_blast} 
        """

# Create SQLite database with all information mikado needs
rule mikado_serialise:
    input:
        mconfig="results/mikado_configure/configuration.yaml",
        blast="results/blast/mikado_prepared.blast.tsv",
        orfs="results/identify_orfs/mikado.orfs.gff3",
        junctions="results/identify_junctions/junctions.bed",
    output:
        db="results/mikado_serialise/mikado.db",
        slog="results/mikado_serialise/serialise.log",
    params:
        extra=config["mikado_serialise_params"]
    conda:
        "../envs/mikado.yaml"
    shell:
        """
        mikado serialise --json-conf {input.mconfig} --orfs {input.orfs} -od results/mikado_serialise/ \
        --junctions {input.junctions} {params.extra}
        """

rule mikado_pick:
    input:
        mconfig="results/mikado_configure/configuration.yaml",
        db="results/mikado_serialise/mikado.db",
    output:
        subloci="results/mikado_pick/mikado.subloci.gff3",
        loci="results/mikado_pick/mikado.loci.gff3"
    params:
        extra=config["mikado_pick_params"]
    conda:
        "../envs/mikado.yaml"
    shell:
        "mikado pick --configuration {input.mconfig} -db {input.db} --subloci_out {output.subloci} -od results/mikado_pick/ {params.extra}"


rule mikado_compare:
    input:
        reference=config["reference_gtf"],
        mikado_out="results/mikado_pick/mikado.loci.gff3"
    output:
        comparison=directory("results/mikado_compare")
    conda:
        "../envs/mikado.yaml"
    shell:
        """
        mikado compare -r {input.reference} --index
        miadko compare -r {input.reference} -p {input.mikado_out} -o results/mikado_compare/compare 
        """
