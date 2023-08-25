# Identify high confidence splice junctions with Portcullis
### SYNTAX CHECKED ###
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
        junctions="results/identify_junctions/junctions.bed",
        get_mikado_input
    output:
        "results/mikado_configure/configuration.yaml"
    params:
        extra=config["mikado_configure_params"]
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
        orfs="results/identify_orfs/mikado.orfs.gff3"
    params:
        extra=config["transdecoder_params"],
    conda:
        "../envs/mikado.yaml"
    log:
        "logs/identify_orfs/TransDecoder.log"
    threads: 1
    shell:
        """
        TransDecoder.LongOrfs -t {input.fasta} --output_dir results/identify_orfs/ {params.extra} 1> {log} 2>&1
        mv transcripts.fasta.transdecoder* results/identify_orfs/
        mv results/identify_orfs/transcripts.fasta.transdecoder.gff3 {output}
        """"


# Run BLAST to get homology data that will help mikado
### SYNTAX CHECKED ###
rule mikado_blast:
    input:
        proteins=config["blast_db"],
        fasta="results/mikado_prepare/mikado_prepared.fasta",
    output:
        prepare_log="results/mikado_blast/blast_prepare.log",
        mikado_blast="results/mikado_blast/mikado_prepared.blast.tsv",
    params:
        extra=config["blastx_params"],
    conda:
        "../envs/mikado.yaml"
    log:
        "logs/mikado_blast/mikado_blast.log"
    threads: 10
    shell:
        """
        makeblastdb -in {input.proteins} -dbtype prot -parse_seqids > {output.prepare_log}
        blastx {params.extra} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qent sstart send evalue bitscore ppos btop" \
        -num_threads {threads} -query {input.fasta} -db {input.proteins} -out {output.mikado_blast} 1> {log} 2>&1
        """

# Create SQLite database with all information mikado needs
### CHECKED SYNTAX ###
rule mikado_serialise:
    input:
        mconfig="results/mikado_configure/configuration.yaml",
        blast="results/blast/mikado_prepared.blast.tsv",
        orfs="results/identify_orfs/mikado.orfs.gff3",
        junctions="results/identify_junctions/junctions.bed",
        blast_db=config["blast_db"]
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
        mikado serialise --json-conf {input.mconfig} --orfs {input.orfs} -od results/mikado_serialise/ \
        --junctions {input.junctions} --xml {input.blast} --blast_targets {input.blast_db} {params.extra} 1> {log} 2>&1
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
        reference=config["reference_gtf"]
    log:
        "logs/make_mikado_list/mikado_list.log"
    threads: 1
    run:

        scallop = params.scallop
        stringtie = params.stringtie
        provided_annotations = params.provided_annotations
        reference = params.reference

        with open(output[0], 'w') as f:

            if scallop["use_scallop"]:

                scallop.pop("use_scallop")

                row = [str(value) for value in scallop.values()]

                row.insert(0, "results/scallop_taco/assembly.gtf")

                f.write('\t'.join(row) + '\n')
            
            if stringtie["use_stringtie"]:

                stringtie.pop("use_stringtie")
                taco = stringtie.pop("use_taco")
                merge = stringtie.pop("use_merge")

                if taco:

                    row = [str(value) for value in stringtie.values()]

                    row.insert(0, "results/stringtie_taco/assembly.gtf")
                
                    f.write('\t'.join(row) + '\n')

                if merge:

                    row = [str(value) for value in stringtie.values()]

                    row.insert(0, "results/stringtie_merge/stringtie_merged.gtf")

                    f.write('\t'.join(row) + '\n')
                

            if ~bool(provided_annotations):
                
                for key, inner_dict in provided_annotations.items():

                    row = [str(value) for value in inner_dict.values()]

                    row.insert(1, key)

                    f.write('\t'.join(row) + '\n')
            
