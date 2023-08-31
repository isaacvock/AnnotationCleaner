# Identify high confidence splice junctions with Portcullis
### SYNTAX CHECKED ###
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
        portcullis prep -t {threads} -o prepare_portcullis {params.extra_prep} {input.fasta} {input.bams} 1> {log} 2>&1
        portcullis junc -t {threads} --orientation {params.orientation} \
        --strandedness {params.strandedness} -o results/identify_junctions/junctions {params.extra_junc} \
        prepare_portcullis/ 1> {log} 2>&1
        """

# Create configuration file
### SYNTAX CHECKED ###
rule mikado_configure:
    input:
        mlist=config["mikado_list"],
        reference=config["genome"],
        proteins=config["blast_db"],
        junctions="results/identify_junctions/junctions.junctions.bed",
        bams=get_mikado_input(),
    output:
        "results/mikado_configure/configuration.yaml"
    params:
        extra=config["mikado_configure_params"],
        scoring=config["mikado_scoring"]
    threads: 4
    conda:
        "../envs/mikado.yaml"
    log:
        "logs/mikado_configure/mikado_configure.log"
    shell:
        """
        mikado configure --list {input.mlist} --scoring {params.scoring} --reference {input.reference} \
        --junctions {input.junctions} -bt {input.proteins} -od results/mikado_configure/ {params.extra} -t {threads} {output} 1> {log} 2>&1
        """  

# Create sorted, non-redundant GTF with all input assemblies
### SYNTAX CHECKED
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
    log:
        "logs/mikado_prepare/mikado_prepare.log"
    threads: 4
    shell:
        "mikado prepare --json-conf {input} -od results/mikado_prepare/ {params.extra} 1> {log} 2>&1"

# Identify open reading frames with TransDecoder
### SYNTAX CHECKED
rule identify_orfs:
    input:
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    output:
        "results/identify_orfs/longest_orfs.gff3"
    params:
        extra=config["transdecoder_params"],
    conda:
        "../envs/transdecoder.yaml"
    log:
        "logs/identify_orfs/TransDecoder.log"
    threads: 1
    shell:
        """
        TransDecoder.LongOrfs -t {input.fasta} --output_dir results/identify_orfs/ {params.extra} 2> {log}
        """

rule predict_orfs:
    input:
        orfs="results/identify_orfs/longest_orfs.gff3",
        transcripts="results/mikado_prepare/mikado_prepared.fasta"
    output:
        "mikado_prepared.fasta.transdecoder.bed"
    params:
        extra=config["transdecoder_predict_params"]
    conda:
        "../envs/transdecoder.yaml"
    log:
        "logs/predict_orfs/TransDecoder.log"
    threads: 1
    shell:
        """
        TransDecoder.Predict -t {input.transcripts} --output_dir results/identify_orfs/ 2> {log}
        """

# Run BLAST to get homology data that will help mikado
### SYNTAX CHECKED ###
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
        makeblastdb -in {input.proteins} -dbtype prot -parse_seqids {params.extra} -out results/mikado_blastdb/mikado_blastdb 1> {log} 2>&1
        """



if config["num_sub"] > 1:
    rule mikado_splitfasta:
        input:
            fasta="results/mikado_prepare/mikado_prepared.fasta",
        output:
            temp(expand("results/mikado_prepare/mikado_prepared.{ID}.fasta", ID = SPLIT_IDS)),
        log:
            "logs/mikado_splitfasta/mikado_splitfasta.log"
        params:
            nsub=config["num_sub"],
        conda:
            "../envs/pyfasta.yaml"
        threads: 1
        shell:
            """
            pyfasta split -n {params.nsub} {input.fasta}
            """

    rule mikado_blastx:
        input:
            db="results/mikado_blastdb/mikado_blastdb.psi",
            fasta="results/mikado_prepare/mikado_prepared.{subID}.fasta",
        output:
            mikado_blast=temp("results/mikado_blast/mikado_prepared.blast.{subID}.tsv"),
        params:
            extra=config["blastx_params"],
        conda:
            "../envs/blast.yaml"
        log:
            "logs/mikado_blastx/mikado_blastx_{subID}.log"
        threads: 6
        shell:
            """
            blastx {params.extra} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
            -num_threads {threads} -query {input.fasta} -db results/mikado_blastdb/mikado_blastdb -out {output.mikado_blast} 2> {log}
            """
            
    rule mikado_mergeblast:
        input:
            expand("results/mikado_blast/mikado_prepared.blast.{ID}.tsv", ID = SPLIT_IDS),
        output:
            "results/mikado_blast/mikado_prepared.blast.tsv",
        conda:
            "../envs/pyfasta.yaml"
        log:
            "logs/mikado_mergeblast/mikado_mergeblast.log"
        threads: 1
        shell:
            """
            awk "NR==1{{print; next}} FNR>1" {input} > {output} 2> {log}
            """

else:

    rule mikado_blastx:
        input:
            db="results/mikado_blastdb/mikado_blastdb.psi",
            fasta="results/mikado_prepare/mikado_prepared.fasta",
        output:
            mikado_blast="results/mikado_blast/mikado_prepared.blast.tsv",
        params:
            extra=config["blastx_params"],
        conda:
            "../envs/blast.yaml"
        log:
            "logs/mikado_blastx/mikado_blastx.log"
        threads: 6
        shell:
            """
            blastx {params.extra} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
            -num_threads {threads} -query {input.fasta} -db results/mikado_blastdb/mikado_blastdb -out {output.mikado_blast} 2> {log}
            """


# Create SQLite database with all information mikado needs
### CHECKED SYNTAX ###
rule mikado_serialise:
    input:
        mconfig="results/mikado_configure/configuration.yaml",
        blast="results/mikado_blast/mikado_prepared.blast.tsv",
        orfs="mikado_prepared.fasta.transdecoder.bed",
        junctions="results/identify_junctions/junctions.junctions.bed",
        blast_db=config["blast_db"],
        transcripts="results/mikado_prepare/mikado_prepared.fasta"
    output:
        db="results/mikado_serialise/mikado.db",
        slog="results/mikado_serialise/serialise.log",
    params:
        extra=config["mikado_serialise_params"]
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

# Identify probable transcripts
### CHECKED SYNTAX ###
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
    log:
        "logs/mikado_serialise/mikado_pick.log"
    threads: 4
    shell:
        "mikado pick --configuration {input.mconfig} -db {input.db} --loci_out mikado.loci.gff3 --subloci_out mikado.subloci.gff3 -od results/mikado_pick/ {params.extra} 1> {log} 2>&1"

# Compare Mikado's to a reference
### Currently cut off from pipeline by target rule ###
rule mikado_compare:
    input:
        reference=config["reference_gtf"],
        mikado_out="results/mikado_pick/mikado.loci.gff3"
    output:
        dummy="results/ignorethisdirectory_mikado/success.txt"
    conda:
        "../envs/mikado.yaml"
    log:
        "logs/mikado_compare/mikado_compare.log"
    threads: 4
    shell:
        """
        mikado compare -r {input.reference} --index 1> {log} 2>&1
        mikado compare -r {input.reference} -p {input.mikado_out} -o results/mikado_compare/compare 1> {log} 2>&1
        touch {output.dummy}
        """



### Create list for mikado (sticking down here to avoid cluttering other rules)
rule make_mikado_list:
    output:
        config["mikado_list"]
    params:
        scallop=config["scallop"],
        stringtie=config["stringtie"],
        provided_annotations=config["provided_annotations"],
    log:
        "logs/make_mikado_list/mikado_list.log"
    threads: 1
    run:

        scallop = params.scallop
        stringtie = params.stringtie
        provided_annotations = params.provided_annotations

        with open(output[0], 'w') as f:

            if scallop["use_scallop"]:

                scallop.pop("use_scallop")
                taco = scallop.pop("use_taco")
                paths = SCALLOP_PATHS

                counter = 1
                for p in paths:
                    
                    row = [str(value) for value in scallop.values()]

                    row.insert(0, str(p))

                    row[1] = row[1] + str(counter)
                    counter = counter + 1

                    f.write('\t'.join(row) + '\n')
            
            if stringtie["use_stringtie"]:

                stringtie.pop("use_stringtie")
                taco = stringtie.pop("use_taco")
                merge = stringtie.pop("use_merge")
                paths = STRINGTIE_PATHS

                counter = 1
                for p in paths:
                    
                    row = [str(value) for value in stringtie.values()]

                    row.insert(0, str(p))

                    row[1] = row[1] + str(counter)
                    counter = counter + 1

                    f.write('\t'.join(row) + '\n')
                

            if ~bool(provided_annotations):
                
                for key, inner_dict in provided_annotations.items():

                    row = [str(value) for value in inner_dict.values()]

                    row.insert(1, key)

                    f.write('\t'.join(row) + '\n')
            
