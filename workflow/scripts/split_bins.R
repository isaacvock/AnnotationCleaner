#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Remove transcripts from annotation that are not supported by read coverage
## relative to intronic background coverage


# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(purrr)
library(furrr)
library(GenomicRanges)
library(optparse)


# Parse command line arguments -------------------------------------------------

args = commandArgs(trailingOnly = TRUE)


option_list <- list(
  make_option(c("-i", "--input", type="character"),
              help = "Path to input flat gtf."),
  make_option(c("-o", "--output", type="chracter"),
              help = "Path to output flat gtf"),
  make_option(c("-s", "--sizelimit"),
              default = 50,
              help = "Largest exon bin size allowed. Larger bins will be split up into smaller bins."),
  make_option(c("-n", "--intronsize"),
              default = 300,
              help = "Largest intron bin size allowed. Larger bins will be split up into smaller bins."),
  make_option(c("-t", "--threads"),
              default = 1,
              help = "Number of threads to use. Binning is parallelized over the intronic/exonic parts using furrr"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


### Function to decrease bin size
#' @param df Input gtf file in data frame format
#' @param size_limit Largest exon bin allowed. Larger bins will be split up into
#' a number of bins as larger or smaller than this.
#' @param debug If TRUE, interactive debugger will be launched
binning <- function(df, size_limit = 200,
                    type = "exonic_part",
                    debug = FALSE){
  
  # Interactive debugger
  if(debug){
    browser()
  }
  
  
  # Only should do this with exonic parts larger than the size limit
  if(df$type == type & df$width > size_limit){
    
    
    # Current bin dimensions
    original_start <- df$start
    original_end <- df$end
    width <- original_end - original_start + 1
    
    ## Create starts and ends for new set of bins
    # Number of bins is such that new bins will be as close to the size limit
    # as possible.
    # Starts should always be the last end + 1, or the original start if the
    # first new bin. Ends should range from the original start + the size limit
    # to the original end. 
    
    num_bins <- ceiling(width/size_limit) + 1
    
    ends <- ceiling(seq(from = original_start + round(size_limit/2),
                        to = original_end,
                        length.out = num_bins))
    
    starts <- c(original_start, ends[1:(length(ends) - 1)] + 1)
    
    # Create final information for gtf
    final_gtf <- data.frame(seqnames = df$seqnames,
                            start = starts,
                            end = ends,
                            strand = df$strand,
                            source = df$source,
                            type = df$type,
                            score = df$phase,
                            gene_id = df$gene_id,
                            transcripts = df$transcripts) %>%
      mutate(width = end - start + 1)
    
    
    return(final_gtf)
    
  }else{
    return(df)
  }
}



# Modify exon bins to increase resolution --------------------------------------

# Load input
flat_gtf_gr <- rtracklayer::import(opt$input)

### Add intronic_parts to gtf

# Split gr into aggregate gene and exonic part
flat_gene_gr <- flat_gtf_gr[mcols(flat_gtf_gr)$type == "aggregate_gene"]
flat_exon_gr <- flat_gtf_gr[mcols(flat_gtf_gr)$type == "exonic_part"]

# Find regions that are different; these are pretty much all introns
difference <- GenomicRanges::setdiff(flat_gene_gr, flat_exon_gr)

# Just to be safe, I will make sure to only consider the differences that
# are also in the aggregate_genes. Everything should overlap but you never know
intron <- GenomicRanges::findOverlaps(difference, flat_gene_gr)

# Construct intronic_part GenomicRanges object
intron_gr <- GRanges(seqnames = seqnames(difference)[queryHits(intron)],
                     ranges = ranges(difference)[queryHits(intron)],
                     strand = strand(difference)[queryHits(intron)])

mcols(intron_gr) <- mcols(flat_gene_gr)[subjectHits(intron),]

mcols(intron_gr)$type <- "intronic_part"

final_gtf_gr <- c(flat_gtf_gr, intron_gr)


### Split up intronic parts into smaller chunks

# Size limit
size_limit <- opt$sizelimit


# Cast to data frame
flat_gtf <- as_tibble(final_gtf_gr)



# Not all parts of gtf need to be worked on
# Set aside aggregate genes or small exon bins
# Will limit number of loop iterations
final_gtf <- flat_gtf %>%
  filter(type == "aggregate_gene" | width <= size_limit)

# Large exonic bins that need to be split up
large_exons <- flat_gtf %>%
  filter(type == "exonic_part" & width > size_limit)

# Large exonic bins that need to be split up
large_introns <- flat_gtf %>%
  filter(type == "intronic_part" & width > size_limit)


### Bin introns and exons

# Set up parallelization
plan(multisession, workers = opt$threads)

# Bin intronic regions
list_of_higher_res_introns <- future_map(seq_len(nrow(large_introns)), ~binning(large_introns[.x, ], 
                                                                                type = "intronic_part",
                                                                                size_limit = opt$sizelimit))
higher_res_introns <- bind_rows(list_of_higher_res_introns)


# Bin exonic regions
list_of_higher_res_exons <- future_map(seq_len(nrow(large_exons)), ~binning(large_exons[.x, ], 
                                                                                     size_limit = opt$intronsize))
higher_res_exons <- bind_rows(list_of_higher_res_exons)


# Combine
higher_res <- bind_rows(higher_res_exons, higher_res_introns)



# Final GTF in data frame form
final_gtf <- bind_rows(final_gtf %>%
                         dplyr::select(-exonic_part_number), 
                       higher_res)

### Prepare for exporting

# Add back exonic part number and create an intronic part number
final_gtf <- final_gtf %>%
  arrange(gene_id, start) %>%
  group_by(gene_id, type) %>%
  mutate(exonic_part_number = ifelse(type != "exonic_part", NA,
                                     1:n()),
         intronic_part_number = ifelse(type != "intronic_part", NA,
                                       1:n()))

# Determine number of digits for leading 0 padding
max_bin_count_e <- max(final_gtf$exonic_part_number[final_gtf$type == "exonic_part"])
num_digits_e <- floor(log10(max_bin_count_e)) + 1

max_bin_count_i <- max(final_gtf$intronic_part_number[final_gtf$type == "intronic_part"])
num_digits_i <- floor(log10(max_bin_count_i)) + 1



# Convert to string and pad with leading 0s
to_string_with_zeros <- function(x, total_chars) {
  sprintf('%0*d', total_chars, x)
}

final_gtf <- final_gtf %>%
  mutate(exonic_part_number = to_string_with_zeros(exonic_part_number, num_digits_e),
         intronic_part_number = to_string_with_zeros(intronic_part_number, num_digits_i))

# Add exon and intron IDs
final_gtf <- final_gtf %>%
  mutate(exon_id = ifelse(type == "exonic_part", 
                          paste0("E", gene_id, exonic_part_number),
                          NA),
         intron_id = ifelse(type == "intronic_part",
                            paste0("I", gene_id, intronic_part_number),
                            NA))

# Export as gtf
final_gr <- GRanges(seqnames = Rle(final_gtf$seqnames),
                    ranges = IRanges(final_gtf$start, end = final_gtf$end),
                    strand = Rle(final_gtf$strand))

mcols(final_gr) <- final_gtf %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)

rtracklayer::export(final_gr,
                    con = opt$output)