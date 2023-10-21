#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Remove unstranded entries of gtf


# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)
library(optparse)


# Parse arguments --------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)


option_list <- list(
  make_option(c("-i", "--input", type="character"),
              help = "Path to input gtf file to be cleaned"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to gtf with unstranded entries removed"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Filter out unstranded --------------------------------------------------------

gtf <- as_tibble(rtracklayer::import(opt$input))

gtf <- gtf %>%
    filter(strand != "*")

# Convert to GenomicRanges object and export
new_GR <- GRanges(seqnames = Rle(gtf$seqnames),
                    ranges = IRanges(gtf$start, end = gtf$end, 
                                    names = 1:nrow(gtf)),
                    strand = Rle(gtf$strand))

mcols(new_GR) <- gtf %>%
    dplyr::select(-seqnames, -start, -end, -strand)

rtracklayer::export(new_GR,
                    con = opt$output)