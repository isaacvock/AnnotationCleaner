import os
import sys
import glob
import itertools
import yaml
from Mikado.utilities import path_join
import Mikado.configuration.configurator
from Mikado.configuration.daijin_configuration import DaijinConfiguration
from Mikado.configuration.configuration import MikadoConfiguration
import subprocess
import gzip
from snakemake import logger as snake_logger
import pkg_resources
import Mikado
import dataclasses
import marshmallow_dataclass

prod_conda = pkg_resources.resource_filename("Mikado.daijin", os.path.join("envs", "prodigal.yaml"))
swissprot = "uniprot_sprot_plants.fasta.gz"
swissprot_noat = "uniprot_sprot_plants.not_at.fasta"

DBs=[swissprot]
zipDBs=["{0}.gz".format(db) for db in DBs]

configname = "configuration.yaml"
config_command = "mikado configure --use-transdecoder --seed 1149199520 --list list.txt --reference chr5.fas.gz --mode permissive --daijin -t 3 \
        --subloci-out mikado.subloci.gff3 --scoring plant.yaml --junctions junctions.bed -bt {swiss} -bc 1 {configname}".format(configname=configname, swiss=swissprot)
toml_name = "configuration.toml"
toml_command_template = "mikado configure --use-transdecoder --seed 1149199520 --list list.txt --reference chr5.fas.gz --mode permissive {daijin} -t 3 \
        -od {outdir}  --subloci-out mikado.subloci.gff3 --scoring plant.yaml --junctions junctions.bed -bt {swiss} -bc 1 {configname}"

toml_command = toml_command_template.format(configname=toml_name, swiss=swissprot, daijin="--daijin", outdir="Daijin")

if not os.path.exists(configname):
    snake_logger.info("Creating the configuration file")
    snake_logger.info(config_command)
    subprocess.call(config_command, shell=True)

snake_logger.info("Reading the configuration file")
dschema = MikadoConfiguration.Schema()
_config = dschema.dump(Mikado.configuration.configurator.load_and_validate_config(configname))

try:
    config = _config
except Exception as exc:
    snake_logger.error("Error in reading the configuration file: {}\n{}".format(exc, _config))
    os.remove(configname)
    raise

# configfile: "configuration.yaml"

rule complete:    
    input: "compare.stats", "compare_subloci.stats", "compare_input.stats", "check.ok", "manual.ok",
           "daijin_test/mikado.yaml", "check_metrics.ok", "g11.ok", "refmap_check.ok", "refmap_check_pc.ok",
           "external.ok", "alias_check.ok"
    output: touch("finished.ok")


rule complete_no_assemble:
    input: "compare.stats", "compare_subloci.stats", "compare_input.stats", "check.ok", "manual.ok",
           "check_metrics.ok", "g11.ok", "refmap_check.ok", "refmap_check_pc.ok",
           "external.ok", "alias_check.ok"
    output: touch("finished_no_daijin.ok")


rule daijin_configure:
    input:
      prots=swissprot,
      genome="chr5.fas.gz"
    output: "daijin.toml"
    shell: "daijin configure -as stringtie scallop -lal gmap --scheduler local -al \
hisat gsnap --sample-sheet samples.txt -o {output} -g {input.genome} -od daijin_test \
 --prot-db {input.prots} --scoring plant.yaml"


rule daijin_assemble:
    input:
        conf=rules.daijin_configure.output
    output:
        conf="daijin_test/mikado.yaml"
    threads: 4
    message: "daijin assemble -nd --nolock --threads 2 --cores 4 --jobs 2 daijin.toml"
    shell: "daijin assemble -nd --nolock --threads 2 --cores 4 --jobs 2 {input.conf}"

rule test_json:
    input: db=swissprot, config=configname
    output: touch("{}.ok".format(configname))
    message: "Checking the configuration"
    run:
        try:
            __= Mikado.configuration.configurator.load_and_validate_config(configname)
        except:
            os.remove(configname)
            raise
        # subprocess.call("gunzip -c chr5.fas.gz > chr5.fas", shell=True)
        # shell("touch {output}")

rule daijin_mikado_configure:
    input: rules.test_json.output
    output: toml_name
    shell: toml_command

rule daijin:
    input: "class.gtf", "cufflinks.gtf", "stringtie.gtf", "trinity.gff3", "mikado.bed", rules.daijin_mikado_configure.output, swissprot
    output:
        loci=os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.gff3"),
        sub=os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado.subloci.gff3"),
	    submetrics=os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado.subloci.metrics.tsv"),
        db=os.path.join("Daijin", "5-mikado", "mikado.db"),
        # mono=os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado.monoloci.gff3"),
        prep=os.path.join("Daijin", "5-mikado", "mikado_prepared.gtf"),
        prep_fasta=os.path.join("Daijin", "5-mikado", "mikado_prepared.fasta"),
	prep_fai=os.path.join("Daijin", "5-mikado", "mikado_prepared.fasta.fai"),
        perm_stats=os.path.join("Daijin", "5-mikado", "pick", "permissive", "mikado-permissive.loci.stats"),
        stats=os.path.join("Daijin", "5-mikado", "pick", "comparison.stats"),
        marker=os.path.join("Daijin", "5-mikado", "all.done")
    message: "daijin mikado --jobs 2 --cores 4 --threads 2 -nd --nolock configuration.toml"
    # shell: "daijin mikado --jobs 2 --cores 4 --threads 2 -nd --nolock configuration.yaml"
    shell: "daijin mikado --jobs 2 --cores 4 --threads 2 -nd --nolock configuration.toml"


rule mikado_configure_manual:
     input: "class.gtf", "cufflinks.gtf", "stringtie.gtf", "trinity.gff3", "mikado.bed", "configuration.toml", swissprot
     output: "configuration.manual.toml"
     shell: toml_command_template.format(configname="{output}", swiss=swissprot, daijin="", outdir="manual")


rule mikado_prepare_manual:
    input: rules.mikado_configure_manual.output
    output:
      gtf="manual/mikado_prepared.gtf",
      fasta="manual/mikado_prepared.fasta"
    shell: "mikado prepare --configuration {input}"


rule mikado_serialise_manual:
    input:
      gtf=rules.mikado_prepare_manual.output.gtf,
      fasta=rules.mikado_prepare_manual.output.fasta,
      config=rules.mikado_configure_manual.output
    output: "manual/mikado.db"
    shell: "mikado serialise --configuration {input.config}"


rule mikado_pick_manual:
    input:
      config=rules.mikado_configure_manual.output,
      gtf=rules.mikado_prepare_manual.output.gtf,
      db=rules.mikado_serialise_manual.output
    output:
      gff="manual/mikado.loci.gff3",
      ok=touch("manual.ok")
    shell: "mikado pick --configuration {input.config}"


rule index_reference:
    input: reference="reference.gff3"
    output: "reference.gff3.midx"
    log: "index.log"
    message: """mikado compare -r {input[reference]} --index --log {log}"""
    shell: """mikado compare -r {input[reference]} --index --log {log}"""

rule compare:
    input: reference="reference.gff3", prediction=rules.daijin.output.loci, midx=rules.index_reference.output
    output: "compare.stats", "compare.tmap", "compare.refmap"
    log: "compare.log"
    message: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare -l {log}"""
    shell: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare -l {log}"""
    
rule compare_input:
    input: reference="reference.gff3", prediction=rules.daijin.output.prep, midx=rules.index_reference.output
    output:
        stats="compare_input.stats",
        tmap="compare_input.tmap",
        refmap="compare_input.refmap"
    log: "compare_input.log"
    message: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_input -l {log}"""
    shell: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_input -l {log}"""

rule compare_subloci:
    input: reference="reference.gff3", prediction=rules.daijin.output.sub, midx=rules.index_reference.output
    output: "compare_subloci.stats", "compare_subloci.tmap", "compare_subloci.refmap"
    log: "compare_subloci.log"
    message: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_subloci -l {log}"""
    shell: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_subloci -l {log}"""

rule compare_subloci_pc:
    input: reference="reference.gff3", prediction=rules.daijin.output.sub, midx=rules.index_reference.output
    output:
        stats="compare_subloci_pc.stats",
        tmap="compare_subloci_pc.tmap",
        refmap="compare_subloci_pc.refmap"
    log: "compare_subloci_pc.log"
    message: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_subloci_pc -l {log}"""
    shell: """mikado compare -r {input[reference]} -p {input[prediction]} -o compare_subloci_pc -l {log}"""

rule check_refmap:
    input:
        refmap=rules.compare_input.output.refmap
    output: touch("refmap_check.ok")
    run:
        import pandas as pd
        refmap = pd.read_csv(input["refmap"], delimiter="\t", index_col=0)
        assert refmap.location.str.contains("^Chr5:", regex=True).all()
        # Account for the non-protein coding
        assert refmap.loc[~refmap.index.str.contains("AT5G66650")].ccode.str.contains(r"^[=|_]$", regex=True).all()
        assert refmap.loc[refmap.index.str.contains("AT5G66650")].ccode.str.contains(r"^f,[=|_]$", regex=True).all()
        assert (refmap.tid.index == refmap.tid.str.replace(r"^at_", "", regex=True)).all(), refmap.tid

rule check_refmap_pc:
    input:
        refmap=rules.compare_subloci_pc.output.refmap
    output: touch("refmap_check_pc.ok")
    run:
        import pandas as pd
        refmap = pd.read_csv(input["refmap"], delimiter="\t", index_col=0)
        assert refmap.location.str.contains(r"^Chr5:", regex=True).all()
        assert refmap.loc[~refmap.index.str.contains("AT5G66650")].ccode.str.contains(r"^[=|_]$", regex=True).all()
        checker = refmap.loc[~refmap.index.str.contains("AT5G66650")].tid.str.replace(r"^at_", "", regex=True)
        assert (checker == checker.index).all()


rule check_logs:
    input:
      pick_out=rules.daijin.output.loci,
      compare_out=rules.compare.output,
      compare_input_out=rules.compare_input.output,
      compare_subloci_out=rules.compare_subloci.output
    output: touch("check.ok")
    run:
      for inp_file in ["Daijin/5-mikado/pick/permissive/mikado-permissive.pick.log",
                       "Daijin/5-mikado/mikado_prepare.log",
                       "Daijin/5-mikado/mikado_serialise.log",
                       "compare.log",
		       "compare_subloci.log", "compare_input.log"]:
         inp_handle = open(inp_file, "rt")
   	 for line in inp_handle:
      	   if "Error" in inp_handle:
      	     raise ValueError(inp_file)
      touch(output[0])

rule check_pick:
    input:
      metrics="Daijin/5-mikado/pick/permissive/mikado.subloci.metrics.tsv",
      fai="Daijin/5-mikado/mikado_prepared.fasta.fai"
    output: "check_metrics.ok"
    run:
      fai = len([_ for _ in open(input["fai"])])
      import pandas as pd
      metrics = pd.read_csv(input["metrics"], delimiter="\t")
      found_split = metrics[metrics.columns[0]].str.contains(r".split").any()
      no_dups = (metrics[metrics.columns[0]].drop_duplicates().shape[0] == metrics.shape[0])
      total = metrics[metrics.columns[0]].str.replace("\.split[0-9]*", "", regex=True).drop_duplicates().shape[0]
      found_excluded = metrics[metrics.columns[2]].str.contains("excluded").any()
      if not found_excluded or not found_split or total != fai or not no_dups:
        raise ValueError("Something went wrong with Mikado pick, please check")
      open(output[0], "wt")


rule test_external_kal_index:
    input:
      fasta=rules.daijin.output.prep_fasta
    output:
      index=os.path.join("Daijin", "5-mikado", "mikado_prepared.idx")
    conda: "kallisto.yaml"
    threads: 1
    shell: "kallisto index -i {output.index} {input.fasta}"

rule test_external_kal_quant:
    input:
      index=rules.test_external_kal_index.output.index,
      r1="ERR588038.R1.fq.gz",
      r2="ERR588038.R2.fq.gz"
    output:
      kal_tsv=os.path.join("Daijin", "5-mikado", "kallisto", "abundance.tsv"),
      data=os.path.join("Daijin", "5-mikado", "data.txt")
    params:
      folder=os.path.join("Daijin", "5-mikado", "kallisto")
    conda: "kallisto.yaml"
    threads: 2
    shell: """kallisto quant -t {threads} -o {params.folder} -i {input.index} {input.r1} {input.r2} && \
cut -f 1,5 {output.kal_tsv} | sed 's/target_id/tid/' > {output.data}"""


rule test_external_kal_serialise:
    input:
      db=rules.daijin.output.db,
      data=rules.test_external_kal_quant.output.data,
    output:
      db=os.path.join("Daijin", "5-mikado", "mikado_external.db"),
      check=touch(os.path.join("Daijin", "5-mikado", "mikado_external.check.ok"))
    params:
      db="mikado_external.db"
    threads: 2
    log: os.path.join("Daijin", "5-mikado", "mikado_serialise_external.log")
    shell: """
    mikado serialise --xml=Daijin/5-mikado/blast/xmls --blast_targets=Daijin/5-mikado/blast/index/blastdb-proteins.fa \
    --start-method=spawn --transcripts=Daijin/5-mikado/mikado_prepared.fasta \
    --genome_fai=Daijin/5-mikado/chr5.fas.gz.fai --json-conf=configuration.yaml \
    --external {input.data} -nsa --force \
    --orfs=Daijin/5-mikado/transdecoder/transcripts.fasta.transdecoder.bed \
    -od Daijin/5-mikado --procs={threads} -l {log} {params.db};
    if [[ $(sqlite3 {output.db} "select count(*) > 0 from external") != 1 ]]; then exit 1; else exit 0; fi"""


rule test_external_kal_pick:
    input:
      db=rules.test_external_kal_serialise.output.db,
      gtf=rules.daijin.output.prep,
      scoring="plant_external.yaml"
    output:
      loci=os.path.join("Daijin", "5-mikado", "pick", "external", "mikado-permissive.loci.gff3"),
      scores=os.path.join("Daijin", "5-mikado", "pick", "external", "mikado-permissive.loci.scores.tsv"),
      subscores=os.path.join("Daijin", "5-mikado", "pick", "external", "mikado.subloci.scores.tsv")
    params:
      outdir=os.path.join("Daijin", "5-mikado", "pick", "external"),
      loci_out="mikado-permissive.loci.gff3"
    log: os.path.join("Daijin", "5-mikado", "pick", "external", "mikado-permissive.log")
    threads: 2
    message: """mikado pick --scoring-file {input.scoring} --source Mikado_permissive \
    --mode=permissive --procs=2 --start-method=spawn \
    --json-conf=configuration.yaml -od {params.outdir} \
     -l {log} --loci-out {params.loci_out} -lv INFO -db {input.db} {input.gtf}"""
    shell: """mikado pick --scoring-file {input.scoring} --source Mikado_permissive \
    --mode=permissive --procs=2 --start-method=spawn \
    --json-conf=configuration.yaml -od {params.outdir} \
     -l {log} --loci-out {params.loci_out} -lv INFO -db {input.db} {input.gtf}"""

rule check_external_pick:
    input:
      scores=rules.test_external_kal_pick.output.scores,
      subscores=rules.test_external_kal_pick.output.subscores,
    output: touch("external.ok")
    run:
        import pandas as pd
        scores = pd.read_csv(input["scores"], delimiter="\t")
        assert "external.tpm" in scores.columns
        assert scores["external.tpm"].max() > 0
        scores = pd.read_csv(input["subscores"], delimiter="\t")
        assert "external.tpm" in scores.columns
        assert scores["external.tpm"].max() > 0

rule check_pick_confusing_alias:
    input:
       db=rules.daijin.output.db,
       prepare="mikado_prepared.conf_alias.gtf"
    output: touch("alias_check.ok")
    log: os.path.join("Daijin", "5-mikado", "pick", "confusing_alias", "pick.log")
    params:
       outdir=os.path.join("Daijin", "5-mikado", "pick", "confusing_alias")
    threads: 2
    shell: "mikado pick --only-reference-update  --source Mikado_permissive --mode=permissive \
    --procs={threads} --start-method=spawn \
    --json-conf=configuration.yaml -od {params.outdir} \
    -l {log} --loci-out mikado-permissive.loci.gff3 -lv INFO \
    -db {input.db} {input.prepare}"

rule test_g11_prodigal:
    input:
      transcripts=rules.daijin.output.prep_fasta
    output: os.path.join("Daijin", "5-mikado", "prodigal", "transcripts.g11.gff3")
    threads: 1
    log: os.path.join("daijin_logs", "prodigal.g11.log")
    conda: prod_conda
    shell: "prodigal -g 11 -i {input.transcripts} -o {output} -f gff 2> {log} > {log}"


rule test_g11_serialise:
    input:
      db=rules.daijin.output.db,
      prodi=rules.test_g11_prodigal.output
    output: os.path.join("Daijin", "5-mikado", "mikado.g11.db")
    threads: 2
    log: os.path.join("Daijin", "5-mikado", "mikado_serialise.g11.log")
    shell: "mikado serialise --xml=Daijin/5-mikado/blast/xmls \
    --blast_targets=Daijin/5-mikado/blast/index/blastdb-proteins.fa --start-method=spawn \
    --transcripts=Daijin/5-mikado/mikado_prepared.fasta --genome_fai=Daijin/5-mikado/chr5.fas.gz.fai \
    --json-conf=configuration.yaml --force --orfs={input.prodi} \
    -od Daijin/5-mikado --procs=2 -l {log} {output}"

rule test_g11_pick:
    input:
      db=rules.test_g11_serialise.output
    output:
      loci=os.path.join("Daijin", "5-mikado", "pick", "g11_permissive", "mikado-permissive.loci.gff3"),
      subloci=os.path.join("Daijin", "5-mikado", "pick", "g11_permissive", "mikado-permissive.subloci.gff3")
    params:
      loci="mikado-permissive.loci.gff3",
      subloci="mikado-permissive.subloci.gff3",
      outdir=os.path.join("Daijin", "5-mikado", "pick", "g11_permissive")
    threads: 2
    log: os.path.join("Daijin", "5-mikado", "pick", "g11_permissive", "mikado-permissive.pick.log")
    shell: "mikado pick --source Mikado_permissive --mode=permissive --procs=2 \
    --start-method=spawn --json-conf=configuration.yaml -od {params.outdir} \
    -l {log} --loci-out {params.loci} --subloci-out {params.subloci} \
    -lv INFO -db {input.db} Daijin/5-mikado/mikado_prepared.gtf"


rule test_g11_gffread:
    input:
      subloci=rules.test_g11_pick.output.subloci,
      genome="chr5.fas.gz"
    output:
      subpep=os.path.join("Daijin", "5-mikado", "pick", "g11_permissive", "mikado-permissive.subloci.pep.fasta"),
      fasta=temp("chr5.fas")
    conda: "gffread.yaml"
    shell: "gunzip -c {input.genome} > {output.fasta} && gffread -S -g {output.fasta} -y {output.subpep} {input.subloci}"


rule check_g11:
    input:
      subpep=rules.test_g11_gffread.output.subpep
    output: touch("g11.ok")
    run:
        with open(input.subpep) as subpep:
            for line in subpep:
                if line[0] != ">" and "*" in line:
                    print(line)
                    raise ValueError(line)


rule clean:
    run:
        import shutil
        for filename in itertools.chain(
                glob.glob("*.ok"), glob.glob("uniprot*.fasta.p*"), glob.glob("daijin_*.stats"),
                glob.glob("*manual*"), glob.glob("*midx"), glob.glob("*fai"), glob.glob("daijin*yaml"),
                glob.glob("uniprot*fasta"), glob.glob("*loci*"),
                ["daijin.yaml", "Daijin", "daijin_logs", "daijin_test"],
                ["mikado_prepared.gtf", "mikado_prepared.fasta", "mikado.prodigal.gff3"],
                glob.glob("compare*"), glob.glob(config["db_settings"]["db"]),
                glob.glob("*.log"), glob.glob("*xml"), ["chr5.fas"],
                ["configuration.yaml", "configuration.toml", "daijin.toml", "daijin.yaml"]):
            if os.path.exists(filename) and filename not in (".", ".."):
                shutil.rmtree(filename) if os.path.isdir(filename) else os.remove(filename)

rule clean_crumbs:
    run:
        for filename in itertools.chain(["finished.ok"], glob.glob("mikado*loci*"),
                                        glob.glob("compare*")):
            if os.path.exists(filename):
                os.remove(filename)
