#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Remove transcripts from annotation that are not supported by read coverage
## relative to intronic background coverage


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
              help = "Largest exon bin size allowed. Larger bins will be split up into smaller bins."))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.



# Modify exon bins to increase resolution --------------------------------------

# Size limit
size_limit <- opt$sizelimit

# Load input
flat_gtf <- as_tibble(rtracklayer::import(opt$input))


### Function to decrease bin size
#' @param df Input gtf file in data frame format
#' @param size_limit Largest exon bin allowed. Larger bins will be split up into
#' a number of bins as larger or smaller than this.
#' @param debug If TRUE, interactive debugger will be launched
exon_binning <- function(df, size_limit = 200,
                         debug = FALSE){
  
  # Interactive debugger
  if(debug){
    browser()
  }
  
  
  # Only should do this with exonic parts larger than the size limit
  if(df$type == "exonic_part" & df$width > size_limit){
    
    
    # Current bin dimensions
    original_start <- df$start
    original_end <- df$end
    width <- original_end - original_start
    
    ## Create starts and ends for new set of bins
    # Number of bins is such that new bins will be as close to the size limit
    # as possible.
    # Starts should always be the last end + 1, or the original start if the
    # first new bin. Ends should range from the original start + the size limit
    # to the original end. 
    
    num_bins <- ceiling(width/size_limit)
    
    ends <- ceiling(seq(from = original_start + size_limit,
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
      mutate(width = end - start)
    
    
    return(final_gtf)
    
  }else{
    return(df)
  }
}


# Not all parts of gtf need to be worked on
# Set aside aggregate genes or small exon bins
# Will limit number of loop iterations
final_gtf <- flat_gtf %>%
  filter(type == "aggregate_gene" | width <= size_limit)

# Large exonic bins that need to be split up
large_exons <- flat_gtf %>%
  filter(type == "exonic_part" & width > size_limit)

# Loop over each overly large bin
for(i in 1:nrow(large_exons)){
  
  if(i == 1){
    
    # Initialize new data frame
    higher_res <- exon_binning(large_exons[i,], size_limit,
                               debug = FALSE)
    
    
  }else{
    
    # Add another row to new data frame
    next_row <- exon_binning(large_exons[i,], size_limit)
    
    higher_res <- bind_rows(higher_res, next_row)
    
  }
  
  
}

# Final GTF in data frame form
final_gtf <- bind_rows(final_gtf %>%
                         dplyr::select(-exonic_part_number), 
                       higher_res)

# Add back exonic part number
final_gtf <- final_gtf %>%
  group_by(gene_id, type) %>%
  mutate(exonic_part_number = ifelse(type == "aggregate_gene", NA,
                                     1:n()))

# Determine number of digits for leading 0 padding
max_bin_count <- max(final_gtf$exonic_part_number[final_gtf$type == "exonic_part"])
num_digits <- floor(log10(max_bin_count)) + 1


# Convert to string and pad with leading 0s
to_string_with_zeros <- function(x, total_chars) {
  sprintf('%0*d', total_chars, x)
}

final_gtf <- final_gtf %>%
  mutate(exonic_part_number = to_string_with_zeros(exonic_part_number, num_digits))


# Export as gtf
final_gr <- GRanges(seqnames = Rle(final_gtf$seqnames),
                    ranges = IRanges(final_gtf$start, end = final_gtf$end),
                    strand = Rle(final_gtf$strand))

mcols(final_gr) <- final_gtf %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)

rtracklayer::export(final_gr,
                    con = opt$output)