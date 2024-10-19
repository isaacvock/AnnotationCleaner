#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Filter out super short transcripts from final merged transcript, which are
## instances of SpliceWiz novel transcripts that mess with StringTie. Other 
## filtering can also be implemented here in the fugure


# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(optparse)


# Parse command line arguments -------------------------------------------------

args = commandArgs(trailingOnly = TRUE)


option_list <- list(
  make_option(c("-i", "--input", type="character"),
              help = "Path to input flat gtf."),
  make_option(c("-o", "--output", type="chracter"),
              help = "Path to output flat gtf"),
  make_option(c("-s", "--sizelimit"),
              default = 200,
              help = "Smallest transcript allowed; all else get filtered"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Filter annotation ------------------------------------------------------------


# Load
gtf <- rtracklayer::import(opt$input)

# Find transcripts of an acceptable length
Transcripts_keep <- as_tibble(gtf) %>%
  filter(type == "exon" & strand != "*") %>%
  group_by(transcript_id) %>%
  summarise(length = sum(width)) %>%
  filter(length >= opt$sizelimit)

# Filter out short transcripts
gtf_filter <- gtf[mcols(gtf)$transcript_id %in% unique(Transcripts_keep$transcript_id)]

# Export
rtracklayer::export(gtf_filter, opt$output)
