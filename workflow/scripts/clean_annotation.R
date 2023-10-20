### PURPOSE OF THIS SCRIPT
## Remove transcripts from annotation that are not supported by read coverage
## relative to intronic background coverage


# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)
library(optparse)


# Parse command line arguments -------------------------------------------------

args = commandArgs(trailingOnly = TRUE)


option_list <- list(
  make_option(c("-r", "--reference", type="character"),
              help = "Path to input gtf file to be cleaned"),
  make_option(c("-f", "--flatref", type = "character"),
              help = "Path to flattened input gtf file"),
  make_option(c("-o", "--output", type="character"),
              help = 'Path to output cleaned gtf file'),
  make_option(c("-e", "--exonic", type="character"),
              default = "",
              help = "Path to exonic counts csv from HTSeq"),
  make_option(c("-b", "--bins", type = "character"),
              default = "",
              help = "Path to exonic bin counts csv from HTSeq"),
  make_option(c("-t", "--total", type = "character"),
              default = "",
              help = "Path to total counts csv from HTSeq"),
  make_option(c("-d", "--directory", type = "character"),
              help = "Path to directory containing all HTSeq count files"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Set opt for testing purposes
opt <- list("G:/Shared drives/Matthew_Simon/IWV/Annotations/Hogg_annotations/stringtie_justin/merged_stringtie_wa.gtf",
         "G:/Shared drives/Matthew_Simon/IWV/Annotations/Hogg_annotations/stringtie_justin/merged_stringtie_wa_flat.gtf",
         "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Testing/cleaned.gtf",
         "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_exonic.csv",
         "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_exonbin.csv",
         "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_total.csv",
         "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/")

names(opt) <- c("reference", "flatref", "output", "exonic", "bins", 
                "total", "directory")

# Source functions for testing purposes
source("C:/Users/isaac/Documents/Simon_Lab/AnnotationPruneR/scripts/Pruners.R")

# Functions called -------------------------------------------------------------


# Prune annotation -------------------------------------------------------------

### Extract Annotations used

# Original
gtf <- as_tibble(rtracklayer::import(opt$reference))

# Flattened + exon IDs
flat_gtf <- as_tibble(rtracklayer::import(opt$flatref))


### Make feature dictionary to relate all feature types to each other
feature_dict <- flat_gtf %>%
  dplyr::filter(type == "exonic_part") %>%
  dplyr::select(gene_id, transcripts, exon_id) %>%
  dplyr::distinct() %>%
  dplyr::mutate(GF = gene_id,
                XF = transcripts,
                all_EF = exon_id) %>%
  dplyr::select(GF, XF, all_EF)

### Extract and process HTSeq output

## Exonic bin quantification
# Due to limitations of HTSeq, this is currently a combination of reads from
# mature and premature RNA

# Get sample IDs
if(opt$bins == ""){
  
  files <- list.files(path = opt$directory, pattern = "exonbin")
  
}else{
  
  files = basename(opt$bins)
  
}


samps <- gsub("_exonbin.csv", "", files)

# Create necessary combined table of all relevant quantifications
# for all samples.
for(i in seq_along(samps)){
  
  ### Exon bin quantification
  
  # File to extract
  filename <- paste0(opt$directory, "/", samps[i], "_exonbin.csv")
  
  # Import file and modify column names
  exonbins_temp <- fread(filename)
  colnames(exonbins_temp) <- c("all_EF", "rname", "reads")
  
  # Filter out last couple meta information rows
  exonbins_temp <- exonbins_temp %>%
    filter(!grepl("__", all_EF))
  
  # Add GF info
  exonbins_temp <- exonbins_temp %>%
    inner_join(feature_dict %>%
                 dplyr::select(GF, all_EF), by  = "all_EF")
  
  # Add mutation rate and sample ID columns
  exonbins_temp <- exonbins_temp %>%
    mutate(mutrate = 0,
           sample = samps[i])
  
  # Add it to growing final table
  if(i == 1){
    
    exonbins <- exonbins_temp
    
  }else{
    
    exonbins <- bind_rows(exonbins, exonbins_temp)
    
  }
  
  ### Exonic quantification
  
  # File to extract
  filename <- paste0(opt$directory, "/", samps[i], "_exonic.csv")
  
  # Import file and modify column names
  exonic_temp <- fread(filename)
  colnames(exonic_temp) <- c("GF", "rname", "exonic_reads")
  
  # Filter out last couple meta information rows
  exonic_temp <- exonic_temp %>%
    filter(!grepl("__", GF)) %>%
    mutate(sample = samps[i])
  
  # Add it to growing final table
  if(i == 1){
    
    exonic <- exonic_temp
    
  }else{
    
    exonic <- bind_rows(exonic, exonic_temp)
    
  }
  
  
  ### Total gene quantification
  
  # File to extract
  filename <- paste0(opt$directory, "/", samps[i], "_total.csv")
  
  # Import file and modify column names
  gene_temp <- fread(filename)
  colnames(gene_temp) <- c("GF", "rname", "total_reads")
  
  # Filter out last couple meta information rows
  gene_temp <- gene_temp %>%
    filter(!grepl("__", GF)) %>%
    mutate(sample = samps[i])
  
  # Add it to growing final table
  if(i == 1){
    
    gene <- gene_temp
    
  }else{
    
    gene <- bind_rows(gene, gene_temp)
    
  }
  
}

### Make dataframes to pass to pruner

# Intronic coverage
intronic_background <- inner_join(gene, exonic, 
                                  by = c("GF", "rname", "sample")) %>%
  mutate(reads = total_reads - exonic_reads) %>%
  filter(reads > 0) %>%
  dplyr::select(GF, rname, reads, sample) %>%
  dplyr::mutate(mutrate = 0)

# Map each sample to an ID corresponding to treatment condition
if(opt$bins == ""){
  
  exp_id = tibble(Exp_ID = c(3, 3, 3, 3, 3, 3,
                             2, 2, 2, 2, 2, 2, 
                             1, 1, 1, 1, 1, 1),
                  sample = unique(intronic_background$sample))

}else{
  
  exp_id = tibble(Exp_ID = 1,
                  sample = unique(intronic_background$sample))
  
  
}

### Score exons
dir <- dirname(opt$output)
score_exons(exonbins, flat_gtf = flat_gtf,
            gtf = gtf, dir = dir,
            exp_ids = exp_id, intronic_background = intronic_background,
            debug = FALSE)

### Prune annotation

EF_to_TF_file <- paste0(dir, "/", "EF_to_TF.csv")
crap_exons_file <- paste0(dir, "/", "unsupported_exons.csv")
exon_check_file <- paste0(dir, "/", "exon_check.csv")

EF_to_TF <- fread(EF_to_TF_file)
crap_exons <- fread(crap_exons_file)
flat_annotation <- flat_gtf
gtf <- gtf
output_path <- opt$output
exon_check <- fread(exon_check_file)

# # Check for which transcripts are garbage
# # Also ignore nuking of first or last exonic bin as this can be 
# # more of a end definition problem than a transcript call problem
# check <- clean_annotation(EF_to_TF = EF_to_TF,
#                           crap_exons = crap_exons,
#                           flat_annotation = flat_annotation,
#                           gtf = gtf,
#                           exon_check = exon_check,
#                           output_path = output_path,
#                           ConsiderEnds = FALSE,
#                           CreateNewGTF = FALSE,
#                           ReturnTranscriptQuality = TRUE)
# 


# Create new annotation
clean_annotation(EF_to_TF = EF_to_TF,
                 crap_exons = crap_exons,
                 flat_annotation = flat_annotation,
                 gtf = gtf,
                 exon_check = exon_check,
                 output_path = output_path,
                 trim = FALSE)

clean_annotation(EF_to_TF = EF_to_TF,
                 crap_exons = crap_exons,
                 flat_annotation = flat_annotation,
                 gtf = gtf,
                 exon_check = exon_check,
                 output_path = paste0(dirname(opt$output), "/cleaned_trimmed.gtf"),
                 trim = TRUE)

exon_check %>% filter(all_EF == "EMSTRG.8907019")

paste0(dirname(opt$output), "/cleaned_trimmed.gtf")
