    rule mikado_serialise:
        input:
            mconfig="results/mikado_configure/configuration.yaml",
            blast=expand("results/mikado_blast/mikado_prepared.blast.{ID}.tsv", ID = SPLIT_IDS),
            orfs="mikado_prepared.fasta.transdecoder.bed",
            junctions="results/identify_junctions/junctions.junctions.bed",
            blast_db=config["blast_db"],
            transcripts="results/mikado_prepare/mikado_prepared.fasta"
        output:
            db="results/mikado_serialise/mikado.db",
            slog="results/mikado_serialise/serialise.log",
        params:
            extra=config["mikado_serialise_params"],
            blast="--tsv={}".format(os.path.join(BLAST_DIR, "tsv")) if len(BLASTX_TARGET) > 0 else "",
        conda:
            "../envs/mikado.yaml"
        log:
            "logs/mikado_serialise/mikado_serialise.log"
        threads: 4
        shell:
            """
            mikado serialise --json-conf {input.mconfig} --transcripts {input.transcripts} --orfs {input.orfs} -od results/mikado_serialise/ \
            --junctions {input.junctions} --tsv {input.blast} --blast_targets {input.blast_db} {params.extra} 1> {log} 2>&1
            """