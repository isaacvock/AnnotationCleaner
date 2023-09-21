rule mikado_blastdb:
    input:
        proteins=config["blast_db"],
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    output:
        "results/mikado_blastdb/mikado_blastdb.psi"
    params:
        extra=config["makeblastdb_params"],
    conda:
        "../envs/blast.yaml"
    log:
        "logs/mikado_blastdb/mikado_blast.log"
    threads: 1
    shell:
        """
        diamond makedb --in {input.proteins} -d prot -parse_seqids {params.extra} -out results/mikado_blastdb/mikado_blastdb 1> {log} 2>&1
        """    
    
    rule mikado_blastx:
        input:
            db="results/mikado_blastdb/mikado_blastdb.psi",
            fasta="results/mikado_prepare/mikado_prepared.{subID}.fasta",
        output:
            mikado_blast="results/mikado_blast/mikado_prepared.blast.{subID}.tsv",
        params:
            extra=config["blastx_params"],
        conda:
            "../envs/blast.yaml"
        log:
            "logs/mikado_blastx/mikado_blastx_{subID}.log"
        threads: 4
        shell:
            """
            blastx {params.extra} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
            -num_threads {threads} -query {input.fasta} -db results/mikado_blastdb/mikado_blastdb -out {output.mikado_blast} 2> {log}
            """