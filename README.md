# AnnotationCleaner

## Introduction
AnnotationCleaner is designed to automate the construction of ab initio assembled transcriptomes using several different tools. This will allow you to test the impact of assembler choice on downstream applications. Currently, the assemblers implemented are:

1. [Scallop](https://github.com/Kingsford-Group/scallop) + [TACO](https://tacorna.github.io/)
2. [Stringtie + Stringtie-merge](https://github.com/gpertea/stringtie)
3. Stringtie + TACO

Other tools that may be implemented in the future are:

1. [Mikado](https://github.com/EI-CoreBioinformatics/mikado)
2. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) + [Cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/)
3. Cufflinks + TACO

Note, TACO is a tool for merging runs of Scallop, Stringtie, or simliar single-sample assemblers on multiple samples. Mikado is a tool combines transcriptomes from multiple different assemblers to infer a consensus transcriptome.

## Quickstart

Below are abriged instructions for running AnnotationCleaner, specifically for those especially comfortable with adopting new bioinformatic tools are who have experience with the key dependencies of AnnotationCleaner (i.e., the assembler tools listed above and Snakemake). The remaining of the README goes into greater detail about each step of this process.

Steps to run AnnotationCleaner:

0. Acquire bam files and a reference annotation (gtf file). These are the input to Annotation Cleaner
1. Install [Git](https://git-scm.com/downloads), [Snakemake](https://snakemake.readthedocs.io/en/stable/), and [Mamba](https://mamba.readthedocs.io/en/latest/) (or [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)).
2. Install Snakedeploy and deploy AnnotationCleaner (for example, follow instructions for deploying workflows on the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe))
   * Create a directory that you want to run AnnotationCleaner in, then run `snakedeploy deploy-workflow https://github.com/isaacvock/AnnotationCleaner.git . --branch main` from inside that directory.
3. Edit config.yaml file, providing paths to bam files and reference annotation, as well as any optional parameters to the assembly tools used.
4. Run with a call to Snakemake that might look something like: `snakemake --cores all --use-conda`. `--use-conda` will automate installation of all other AnnotationCleaner dependencies, and is thus highly recommended. Otherwise, all dependencies will have to be installed manually and made available for AnnotationCleaner to call.

## Running 

### Step 1: Install basic dependencies

AnnotationCleaner uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies (e.g., HTSeq) are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

Finally, installation of all other dependencies is automated with [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)/[Mamba](https://mamba.readthedocs.io/en/latest/). Conda is a package/environment management system. Mamba is a new, faster, C++ reimplementation of conda. While often associated with Python package management, lots of software, including all of the TimeLapse pipeline dependencies, can be installed with these package managers. They have pretty much the same syntax and can do the same things, so **I highly suggest using Mamba in place of Conda whenever possible**. 

One way to install Mamba is to first install Conda following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can call:

``` bash
conda install -n base -c conda-forge mamba
```
to install Mamba.

A second strategy would be to install Mambaforge, which is similar to something called Miniconda but uses Mamba instead of Conda. I will reproduce the instructions to install Mambaforge below, as this is probably the easiest way to get started with the necessary installation of Mamba. These instructions come from the [Snakemake Getting Started tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html), so go to that link if you'd like to see the full original details:

* For Linux users with a 64-bit system, run these two lines of code from the terminal:

``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
* For Mac users with x86_64 architecture: 
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```
* And for Mac users with ARM/M1 architecture:
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
bash Mambaforge-MacOSX-arm64.sh
```

When asked this question:
``` bash
Do you wish the installer to preprend the install location to PATH ...? [yes|no]
```
answer with `yes`. Prepending to PATH means that after closing your current terminal and opening a new one, you can call the `mamba` (or `conda`) command to install software packages and create isolated environments. We'll be using this in the next step.

AnnotationCleaner requires bam files as input. These bam files can be obtained via alignment of fastq files with any of a number of aligners. [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://github.com/DaehwanKimLab/hisat2) are commonly suggested for this task.

### Step 2: Deploy workflow

This workflow can easily be deployed for your particular usecase with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html). To get started with Snakedeploy, you first need to create a simple conda environment with Snakemake and Snakedeploy:


``` bash
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy
```

Next, create a directory that you want to run bam2bakR in (I'll refer to it as `workdir`) and move into it:
``` bash
mkdir workdir
cd workdir
```

Now, activate the `deploy_snakemake` environment and deploy the workflow as follows:

``` bash
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/AnnotationCleaner.git . --branch main
```

`snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR.git` copies the content of the `config` directory in the bam2bakR Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of the bam2bakR Github repo.

### Step 3: Update config

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  Rep_1: data/bam/WT_replicate_1.bam
  Rep_2: data/bam/WT_replicate_2.bam
  Rep_3: data/bam/WT_replicate_3.bam
```
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `Rep_1`, `Rep_2`, and `Rep_3`. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path is NOT an absolute path, it is relative to the directory that you deployed to (i.e., `workdir` in this example). Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the bam2bakR directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the paths would look something like this:

``` yaml
samples:
  WT_1: ../data/WT_replicate_1.bam
  WT_2: ../data/WT_replicate_2.bam
  WT_ctl: ../data/WT_nos4U.bam
  KO_1: ../data/KO_replicate_1.bam
  KO_2: ../data/KO_replicate_2.bam
  KO_ctl: ../data/KO_nos4U.bam
```
where `../` means navigate up one directory. 

The second parameter you will have to edit is `reference_gtf`, which is an annotation file that will be used to map the assembled transcripts to previously annotated transcripts and genes. The gtf should be from the same genome that you used for alignment of fastq files. 

The rest of the parameters in the config decide which assemblers will be used, and tunes their behavior. For example, as shown in the example config, you will want to specify the strandedness of your library in the `scallop_params` and `stringtie_params` parameters. 

For Scallop, the options of what you can specify for the strandedness (`--library_type`)are `unstranded`, `first`, and `second`. `unstranded` means that the library was not stranded. `first` means that the first read in a read pair is the reverse complement of the original RNA sequence (and the second read is the original RNA sequence). `second` means that the first read in a read pair is the original RNA sequence.

For Stringtie, the options of what you can specify for the strandedness are `--rf` (equivalent of `library_type first` in Scallop) and `fr` (equivalent of `library_type second` in Scallop).

See the documentation for relevant tools to get information about all other parameters that can be set.

### Step 4: Run

Once steps 1-3 are complete, AnnotationCleaner can be run from the directory you deployed the workflow to as follows:

``` bash
snakemake --cores all --use-conda
```
There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I will point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.


## Output

In the `results` directory, you will find the following sub-directories generated by AnnotationCleaner (if you ran both Stringtie and Scallop; only the directories related to the relevant tool will be present if you ran a single one of these):

1. `sorted`: Sorted bam files
2. `separate_scallops`: Scallop gtfs created for each individual sample passed in
3. `scallop_taco`: Output of `taco_run` using all Scallop gtfs as input
4. `scallop_taco_refcomp`: Output of `taco_refcomp` using assembly.gtf present in `scallop_taco` and the provided reference gtf
5. `separate_stringties`: Stringtie gtfs created for each individual sample passed in
6. `stringtie_taco`: Output of `taco_run` using all Stringtie gtfs as input
7. `stringtie_merge`: Output of `stringtie --merge` using all Stringtie gtfs as input
8. `stringtie_taco_refcomp`: Output of `taco_refcomp` using assembly.gtf present in `stringtie_taco` and the provided reference gtf
9. `ignorethisdirectory`: Directory that had to get created as a cheap hack for Snakemake to run TACO with Scallop output (problem is Snakemake needs some output for every rule, which in TACO's case should be a directory, but specifying this output to Snakemake causes Snakemake to create the directory prior to TACO running, which breaks TACO). Contains a singular empty text file
10. `ignorethisdirectory_stringtie`: Directory that had to get created as a cheap hack for Snakemake to run TACO with Stringtie output. Contains a singular empty text file
