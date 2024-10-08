#!/usr/bin/env Rscript
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
  make_option(c("-b", "--bins", type = "character"),
              default = "",
              help = "Path to exonic bin counts from FeatureCounts"),
  make_option(c("-u", "--intronbins", type = "character"),
              default = "",
              help = "Path to intron bin counts from FeatureCounts"),
  make_option(c("-d", "--directory", type = "character"),
              help = "Path to directory containing all FeatureCounts count files"),
  make_option(c("-c", "--fdr", type = "double"),
              default = 0.05,
              help = "Conservativeness of exon bin removal (FDR for true exon bins)"),
  make_option(c("-l", "--floor", type = "double"),
              default = 2.0,
              help = "Minimum factor difference between exonic and intronic RPK"),
  make_option(c("-a", "--readlength", type = "integer"),
              default = 100,
              help = "Total effective read length"),
  make_option(c("-p", "--priorvar", type = "double"),
              default = 2,
              help = "Variance of prior on intronic RPK"),
  make_option(c("-s", "--slope", type = "double"),
              default = 5,
              help = "Non-linear slope parameter for coverage vs. dispersion trend"),
  make_option(c("-i", "--intercept", type = "double"),
              default = 0.01,
              help = "Non-linear intercept parameter for coverage vs. dispersion trend"),
  make_option(c("-y", "--ifactor", type = "double"),
              default = 1,
              help = "Intronic RPK estimate = pmin(median(RPK) + ifactor*mad(RPK), max(RPK))"),
  make_option(c("-g", "--groups", type = "string"),
              default = "",
              help = "Path to table of sample to group ID mappings. If not provided, all samples
              assumed to be from the same group.")  )

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# # Set opt for testing purposes
# opt <- list("G:/Shared drives/Matthew_Simon/IWV/Annotations/Hogg_annotations/stringtie_justin/merged_stringtie_wa.gtf",
#          "G:/Shared drives/Matthew_Simon/IWV/Annotations/Hogg_annotations/stringtie_justin/merged_stringtie_wa_flat.gtf",
#          "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Testing/cleaned.gtf",
#          "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_exonic.csv",
#          "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_exonbin.csv",
#          "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/NMD11j_ctl_1_total.csv",
#          "G:/Shared drives/Matthew_Simon/IWV/FlatStacks/with_htseq/st_justin_v1/quantify/",
#          0.05, 3, 100, 2, 5, 0.01, 10, 50, TRUE, TRUE, TRUE, 10)
# 
# names(opt) <- c("reference", "flatref", "output", "exonic", "bins", 
#                 "total", "directory", "fdr", "floor", "readlength",
#                 "priorvar", "slope", "intercept", "meancutoff", "maxcutoff",
#                 "discardgenes", "ignoreends", "notrimming", "minreads")
# 
# # Source functions for testing purposes
# source("C:/Users/isaac/Documents/Simon_Lab/AnnotationPruneR/scripts/Pruners.R")

# Functions called -------------------------------------------------------------

# Helper functions that I will use on multiple occasions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))


### Function to score exons as likely intronic or not
score_exons <- function(cB, flat_gtf, gtf, dir,
                        exp_ids, intronic_background = NULL,
                        FDR = 0.05, coverage_floor = 2,
                        read_length = 100, tau2 = 2, ifactor=1.5,
                        a1 = 5, a0 = 0.01, sampleID = "",
                        debug = FALSE){
  
  if(debug){
    browser()
  }
  
  
  ### Check if a cB or HTSeq output was provided as input
  cB_cols <- c("sample", "EF", "GF", "XF", "TC", "nT")
  htseq_cols <- c("sample", "all_EF", "GF", "mutrate")
  
  if(all(cB_cols %in% colnames(cB))){
    
    if(!is.null(intronic_background)){
      warning("Both cB and intronic background data frames were provided.
              The intronic background data frame will be overwridden by
              creation of said data frame from the cB data frame.")
    }
    
    input_cB <- TRUE
    
  }else if(all(htseq_cols %in% colnames(cB)) & !is.null(intronic_background)){
    
    input_cB <- FALSE
    
  }else{
    
    stop("!!Input is neither a cB data frame nor properly formatted HTSeq output!!")
    
  }
  
  
  ### Preprocess cB if that is what got passed in
  if(input_cB){
    
    # Filter out garbage and spike-ins (though I think latter shouldn't even be present given genome I used)
    cB <- cB[!grepl("__", GF) & !grepl("dm", rname)]
    
    
    ### Calculate intronic background
    
    # Reads and mutation rates in definitely premature RNA reads
    cB_i <- cB[grepl("__", XF)]
    
    # Coverage and mutational content over introns
    intronic_background <- cB_i[,.(reads = sum(n), mutrate = sum(TC*n)/sum(nT*n)),
                                by = .(sample, GF)]
    
    ### Assign each read to all of its exon bins
    
    # Remove extra text
    cB$EF <- gsub('__ambiguous', '', cB$EF)
    cB$EF <- gsub('\\[', '', cB$EF)
    cB$EF <- gsub('\\]', '', cB$EF)
    
    
    # Map ambiguous to individual bins
    EF_dict <- cB[,c("EF", "EF")] %>%
      dplyr::distinct() 
    
    colnames(EF_dict) <- c("XF", "EF")
    
    EF_dict <- EF_dict %>%
      separate_rows(EF, sep = "\\+") 
    
    colnames(EF_dict) <- c("EF", "all_EF")
    
    
    # Merge
    cB <- cB[, .(n = sum(n)), by = .(sample, EF, GF, XF, TC, nT)]
    
    cB <- inner_join(cB, EF_dict, by = "EF", 
                     relationship = "many-to-many")
    
    cB <- cB[,.(reads = sum(n), mutrate = sum(TC*n)/sum(nT*n)), 
             by = .(sample, all_EF, GF)]
    
    gc(cB_i)
    
    # Get intron lengths to normalize
    intron_sizes <- flat_gtf %>%
      group_by(gene_id) %>%
      summarise(total_length = width[type == "aggregate_gene"],
                exon_length = sum(width[type == "exonic_part"])) %>%
      mutate(intron_length = total_length - exon_length) %>%
      filter(intron_length >= 0)
    colnames(intron_sizes) <- c("GF", "total_length", "exon_length", "intron_length")
    
    intronic_background <- inner_join(intronic_background, intron_sizes, by = c("GF"))
    
    
    # # Let's just filter out 0 intron length stuff 
    # # with the caveat that we can't say anything about validity of its annotation
    #   # ACTUALLY, with hierarchical modeling of intronic RPK estimates, we can!
    # intronic_background <- intronic_background[intron_length > 0]
    
    
    # RPM normalize
    intronic_background <- intronic_background %>%
      mutate(RPK = reads/((intron_length + ((read_length - 1)))/1000))
    
    
    
  }else{
    
    # Get intron lengths to normalize
    intron_sizes <- flat_gtf %>%
      filter(type == "intronic_part") %>%
      dplyr::mutate(intron_length = width,
                    GF = gene_id,
                    all_IF = intron_id) %>%
      dplyr::select(GF, all_IF, intron_length)
    
    # Intronless genes  
    intronless_genes <- flat_gtf %>%
      group_by(gene_id) %>%
      summarise(total_length = width[type == "aggregate_gene"],
                exon_length = sum(width[type == "exonic_part"])) %>%
      mutate(intron_length = total_length - exon_length) %>%
      filter(intron_length == 0) %>%
      mutate(RPK = 0,
             mutrate = 0,
             sample = sampleID) %>%
      dplyr::select(sample, gene_id, RPK, mutrate, intron_length)
      
    colnames(intronless_genes) <- c("sample", "GF", "RPK", "mutrate", "intron_length")
    
    
    
    # RPK normalize
      # If there is more than 1 intronic bin, use the minimum of either the
      # max RPK of all bins or the median RPK plus some constant (user-defined, default of 1.5) times
      # the mad of RPKs (mad calculated on log-scale). Captures uncertainty in the intronic
      # RPK estimate that the negative binomial model downstream of this 
      # does not currently account for.
    intronic_background <- intronic_background %>%
      mutate(RPK = (reads + 1)/((intron_length + ((read_length - 1)))/1000)) %>%
      group_by(sample, GF) %>%
      summarise(RPK = ifelse(n() < 2, 
                             mean(RPK),
                             pmin(median(RPK) + ifactor*stats::mad(RPK), max(RPK))),
                intron_length = sum(intron_length),
                mutrate = mean(mutrate))
    
    
    intronic_background <- bind_rows(intronic_background,
                                     intronless_genes)
    
    
  }
  
  
  
  

  
  ### Now save a map from EFs to transcripts for use later
  
  # Get information mapping exon bin to set of transcripts
  EF_to_TF <- flat_gtf %>%
    dplyr::filter(type == "exonic_part") %>%
    dplyr::select(exon_id, gene_id, transcripts) %>%
    dplyr::distinct() 
  colnames(EF_to_TF) <- c("all_EF", "GF", "XF")
  
  # Map set of transcripts to individual transcripts
  TF_dict <- EF_to_TF %>%
    dplyr::select(XF) %>%
    dplyr::distinct() %>%
    dplyr::mutate(all_XF = XF) %>%
    tidyr::separate_rows(all_XF, sep = "\\+") %>%
    dplyr::distinct()
  
  # Map exon bin to individual transcripts
  EF_to_TF <- inner_join(EF_to_TF, TF_dict, by = "XF", 
                         relationship = "many-to-many")
  
  write_csv(EF_to_TF, 
            file = paste(dir, paste0(sampleID, "EF_to_TF.csv"), sep = "/"))
  write_csv(intronic_background,
            file = paste(dir, paste0(sampleID, "intronic_background.csv"), sep = "/"))
  
  
  ### Calculate coverage over exonic bins
  
  
  
  # Add intronic coverage info
  exon_check <- inner_join(cB, intronic_background %>%
                             dplyr::mutate(intronic_mutrate = mutrate) %>%
                             dplyr::select(sample, GF, intronic_mutrate, RPK, intron_length),
                           by = c("sample", "GF"))
  
  gc(cB)
  gc(EF_to_TF)
  
  ## RPK normalize exonic coverages
  
  
  # Calculate exon lengths
  exon_sizes <- flat_gtf %>%
    filter(!is.na(exon_id) & type == "exonic_part") %>%
    mutate(exon_length = width,
           all_EF = exon_id) %>%
    dplyr::select(all_EF, exon_length)
  
  # Add to exon length to account for read length
  # Rationale is that read of length L can intersect with exon if any of its
  # bases intersect with the exon. In the extreme case, this means 1 base on
  # either end intersecting, and the rest not intersecting. Thus, the effective
  # length = exon length + (read length - 1) 
  # The (read length - 1) is the furthest extent from the 5' end of the exon
  # that the read start can be to still intersect with the exon. 
  #
  # NOTE FOR FUTURE ISAAC (and potentially others): 
  # at first I was fooled into double counting this factor to account for 1 base intersecting
  # with the 5' end of the exon AND the alternate scenario of 1 base intersecting
  # with the 3' end of the exon. This is an inconsistent definition of exon length
  # though. Define the exon length as the number of read start bases that could
  # yield an overlap of the read with the exon. This includes the read starting
  # (read length - 1) bases upstream of the 5' end and concludes at the 3' end of
  # the exon. Accounting for the 3' extent of the read means simultaneously defining
  # exon length as number of start bases that could yield an overlap AND the number
  # of end bases that could yield an overlap. Gotta choose one definition or the other,
  # not both.
  exon_check <- inner_join(exon_check , exon_sizes, by = "all_EF")
  exon_check <- exon_check %>%
    mutate(RPK_exon = reads/((exon_length + ((read_length - 1)))/1000))
  
  ### Estimate exonic vs. intronic coverage trend
  
  # Max exonic coverage and average intronic coverage are used in regression model
  exon_vs_intron <- exon_check %>%
    group_by(GF, sample) %>%
    summarise(intronic_avg = mean(RPK),
              exonic_max = max(RPK_exon))
  
  # Fit linear model and extract parameters of interest
  trend <- exon_vs_intron %>%
    filter(intronic_avg > 0 & exonic_max > 0) %>%
    group_by(sample) %>%
    do(fit = lm(log10(intronic_avg) ~ log10(exonic_max), 
                data = .)) %>%
    rowwise() %>%
    mutate(slope = summary(fit)$coefficients[2,1],
           intercept = summary(fit)$coefficients[1,1]) %>%
    dplyr::select(-fit)
  
  
  
  
  ### Calculate regularized intronic RPK estimate
  
  # Calculate posterior
  regularized <- as_tibble(exon_check) %>%
    inner_join(trend, by = "sample") %>%
    inner_join(exon_vs_intron, by = c("sample", "GF")) %>%
    group_by(sample, GF) %>%
    mutate(mu = RPK*((intron_length + (read_length - 1))/1000),
           sig2_1 = 1/(mu*(log(10)^2)),
           sig2_2 = 1/((log(10)^2)*((a1/mu) + a0)),
           sig2 = sig2_1 + sig2_2,
           mu_post_num = (1/sig2)*log10(RPK) + (slope*log10(exonic_max) + intercept)/(mean(sig2)/tau2),
           mu_post_den = (1/sig2) + (tau2/mean(sig2)),
           mu_post = ifelse(intron_length == 0, 
                            (slope*log10(exonic_max) + intercept),
                            mu_post_num/mu_post_den)) %>%
    ungroup() %>%
    mutate(RPK_post = 10^mu_post) %>%
    mutate(RPK_post = ifelse(RPK_post > RPK, RPK_post, RPK)) %>%
    dplyr::select(-mu_post_den,
                  -mu_post_num,
                  -sig2,
                  -sig2_2,
                  -sig2_1,
                  -mu)
  
  ### Negative binomial model to estimate p-values for intronic null
  
  exon_check <- regularized %>%
    mutate(pval = pnbinom(q = reads,
                          mu = RPK_post*((exon_length + ((read_length - 1)))/1000),
                          size = 1/((a1/(RPK_post*((exon_length + ((read_length - 1)))/1000))) + a0),
                          lower.tail = FALSE)) %>%
    ungroup() %>%
    mutate(padj = stats::p.adjust(pval, method = "BH"))
  
  
  ### Figure out which bins are supported to exist
  
  # Make sure it has good coverage
  exon_check <- exon_check %>%
    mutate(pass = ifelse(padj < FDR & RPK_exon > coverage_floor*RPK_post, 
                         TRUE, FALSE))
  
  
  exon_check <- inner_join(exon_check, exp_ids, by = "sample")
  
  write_csv(exon_check,
            file = paste(dir, paste0(sampleID, "exon_check.csv"), sep = "/"))
  
  
  # Make sure it passes filter in at least one sample
  bonafide_exons <- exon_check %>%
    group_by(all_EF, GF, Exp_ID) %>%
    summarise(pass = ifelse(any(pass), TRUE, FALSE)) %>%
    group_by(all_EF, GF) %>%
    summarise(pass = ifelse(any(pass), TRUE, FALSE))
  # Filtered out some 34,000 exon bins once I properly kept all exonic reads
  #   and used a 3*RPK filter, plus Exp_ID based filtering
  # 30,000 with a 2*RPK filter
  # 25,000 with a 1*RPK filters
  # 56,201 with version 3, FDR < 0.05 and RPK > 3x intronic
  
  crap_exons <- bonafide_exons %>%
    dplyr::filter(!pass) %>%
    dplyr::select(all_EF, GF)
  
  real_exons <- bonafide_exons %>%
    dplyr::filter(pass) %>%
    dplyr::select(all_EF, GF)
  
  
  write_csv(crap_exons,
            file = paste(dir, paste0(sampleID, "unsupported_exons.csv"), sep = "/"))
  write_csv(real_exons,
            file = paste(dir, paste0(sampleID, "supported_exons.csv"), sep = "/"))
  
  
}


### Function to process and create new annotation
clean_annotation <- function(EF_to_TF,
                             crap_exons,
                             flat_annotation,
                             gtf,
                             output_path,
                             dir,
                             exon_check = NULL,
                             debug = FALSE){
  
  
  
  if(debug){
    browser()
  }
  
  
  ##### PRUNE
  
  
  # Find intronless genes.
  # Handled differently when trimming ends and determining if "trash
  transcript_lengths <- gtf %>%
    filter(type == "transcript") %>%
    dplyr::select(transcript_id, width) 
  
  intronless_genes <- gtf %>%
    filter(type == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    summarise(exonic_length = sum(width)) %>%
    inner_join(transcript_lengths, by = "transcript_id") %>%
    mutate(intron_length = width - exonic_length) %>%
    filter(intron_length == 0) %>%
    dplyr::select(transcript_id) %>%
    dplyr::distinct()
  
  
  # Flag exonic bins preceding first sj and proceding last sj
  sj_annotated <- flat_annotation %>%
    filter(type == "exonic_part") %>%
    inner_join(EF_to_TF %>%
                 dplyr::mutate(transcripts = XF) %>%
                 dplyr::select(all_XF, transcripts) %>%
                 dplyr::distinct(),
               by = "transcripts",
               relationship = "many-to-many") %>%
    mutate(exonic_part_number = as.numeric(exonic_part_number)) %>%
    group_by(all_XF) %>%
    mutate(followed_by_sj = ifelse(exonic_part_number != max(exonic_part_number), 
                                   ifelse(lead(start, order_by = exonic_part_number) == (end + 1),
                                          FALSE,
                                          TRUE), FALSE),
           preceded_by_sj = ifelse(exonic_part_number != min(exonic_part_number), 
                                   ifelse(lag(end, order_by = exonic_part_number) == (start - 1),
                                          FALSE,
                                          TRUE), FALSE)) 
  # Technically followed_by_sj will be true whenever an exonic bin is followed by
  # a splice junction. I think I want to account for edge case where I split exonic bins
  # into smaller bins and thus want to just identify any bin that is at risk of getting
  # pruned and is after the first splice junction and before the last
  
  
  # Add exonic bin number and prune crap bins
  EF_test <- EF_to_TF %>% #filter(GF == "SNRNP70") %>%
    inner_join(sj_annotated %>% 
                 ungroup() %>%
                 filter(type == "exonic_part") %>%
                 dplyr::mutate(all_EF = exon_id) %>%
                 dplyr::select(start, end, all_EF, exonic_part_number, all_XF, 
                               followed_by_sj, preceded_by_sj) %>%
                 dplyr::distinct(),
               by = c("all_EF", "all_XF"),
               relationship = "many-to-many")
  
  
  ### TO-DO: NEED TO HAVE AN INTRONLESS EDGE CASE HANDLE WHERE I DO A BETTER JOB
  ### OF NOT PRUNING ANY BINS WHICH WOULD INTRODUCE A NEW SPLICE JUNCTION
  
  ### Can I do it all from EF_test?
  # I think so! And it was stupid easy...
  EF_test2 <- EF_test %>%
    # annotate some bins as likely intronic
    mutate(likely_intronic = all_EF %in% crap_exons$all_EF) %>%
    # annotate which bins are prunable
    group_by(GF, all_XF) %>%
    mutate(first_good = ifelse(
      all(likely_intronic), Inf,
      min(exonic_part_number[!likely_intronic])
    ),
    last_good = ifelse(
      all(likely_intronic), 0,
      max(exonic_part_number[!likely_intronic])
    ),
    last_bin_b4_1st_sj = ifelse(!any(followed_by_sj),
                                Inf,
                                min(exonic_part_number[followed_by_sj])),
    first_bin_after_last_sj = ifelse(!any(preceded_by_sj),
                                     1,
                                     max(exonic_part_number[preceded_by_sj]))) %>%
    ungroup() %>%
    mutate(cannot_prune = (first_good > last_bin_b4_1st_sj) |
    (is.infinite(first_good)) |
    (last_good < first_bin_after_last_sj)) %>%
    ungroup() %>%
    # Annotate which bins should be pruned
    mutate(prune = ((exonic_part_number < first_good) |
                      (exonic_part_number > last_good)) & !cannot_prune) %>% # Don't flag any exonic bins as prunable if pruning introduces splice junctions
    group_by(GF, all_XF) %>%
    mutate(problematic = cannot_prune |
             (any(exonic_part_number[likely_intronic] >= unique(last_bin_b4_1st_sj) & exonic_part_number[likely_intronic] <= unique(first_bin_after_last_sj))) |
             is.infinite(unique(first_good))) # Catch transcripts with bad internal bins or no good bins at all
  
  # Save this useful table with details regaring all exons and their prunability
  write_csv(EF_test, file = paste0(dir, "/exon_details.csv"))
  
  EF_test2 %>%
    group_by(GF, all_XF) %>%
    summarise(problematic = all(problematic)) %>%
    dplyr::group_by(problematic) %>%
    dplyr::count()
  
  # Collapse exon bin graph to allow for cute string grepping
  # trick in next step
  check2 <- EF_test2 %>%
    group_by(GF, all_XF) %>%
    mutate(score = sum(prune)) %>%
    dplyr::filter(!prune) %>%
    summarise(collapsed_graph = paste(exonic_part_number, collapse = ","),
              score = mean(score))
  
  
  ## Function to identify if one string is contained in another  
  # Since transcripts are converted to string representation in last step,
  # this is identical to finding instances where one transcript is 
  # a subset of the exonic parts that make up another transcript.
  is_contained <- function(strings) {
    sapply(strings, function(x) {
      any(sapply(strings, function(y) {
        if(x != y) {
          vect1 <- strsplit(x, ",")
          vect2 <- strsplit(y, ",")
          start <- vect1[1]
          end <- vect1[length(vect1)]
          if(!(start %in% vect2)){
            FALSE
          }else if(!(end %in% vect2)){
            FALSE
          }else{
            substart <- which(vect2 == start)
            subend <- which(vect2 == end)
            subvect2 <- vect2[substart:subend]
            if(identical(x, subvect2)){
              TRUE
            }else{
              FALSE
            }
          }
        } else {
          FALSE
        }
      }))
    })
  }
  
  # Filter out modified transcripts that are now subsets of other transcripts
  check2 <- check2 %>%
    group_by(GF) %>%
    filter(!is_contained(collapsed_graph) | score == 0)
  
  # Keep only the maximum score transcript of a set of identical transcripts
  # If multiple max score transcripts in a set, then just choose one at random
  check2 <- check2 %>%
    ungroup() %>%
    group_by(GF, collapsed_graph) %>%
    filter(score == max(score)) %>%
    summarise(all_XF = all_XF[1],
              GF = GF[1],
              collapsed_graph = collapsed_graph[1],
              score = score[1])
  
  
  ## Now clean out annotation
  
  # 1) Remove garbage transcripts
  check <- EF_test2 %>%
    dplyr::filter(!prune) %>%
    inner_join(check2 %>% dplyr::select(GF, all_XF) %>% dplyr::distinct(),
               by = c("GF", "all_XF")) %>%
    group_by(GF, all_XF) %>%
    summarise(new_start_base = start[exonic_part_number == min(exonic_part_number)], # DO I NEED TO ASSESS STRANDEDNESS?
              new_end_base = end[exonic_part_number == max(exonic_part_number)],
              problematic = any(problematic))
  
  
  ### NOTE: CURRENTLY THERE IS A WARNING WHEN I INNER JOIN BECAUSE CHECK HAS
  ### ONE ROW FOR EACH RETAINED EXONIC BIN
  
  # Filter out bogus transcripts and prune ends of real transcripts
  transcript_annotation <- gtf %>%
    filter(type == "transcript") %>%
    inner_join(check %>%
                 ungroup() %>%
                 mutate(transcript_id = all_XF) %>%
                 dplyr::select(transcript_id, new_start_base, new_end_base, problematic),
               by = "transcript_id") %>%
    mutate(start = new_start_base,
           end = new_end_base) %>%
    mutate(width = end - start)  %>%
    dplyr::select(-new_start_base, -new_end_base)
  
  
  ### TO-DO: fix fact that all exons are getting the new start and ends
  # Filter out bogus exons and prune ends of real exons
  exon_annotation <- gtf %>%
    filter(type == "exon") %>%
    inner_join(check %>%
                 ungroup() %>%
                 mutate(transcript_id = all_XF) %>%
                 dplyr::select(transcript_id, new_start_base, new_end_base, problematic),
               by = "transcript_id") %>%
    mutate(start = ifelse(start < new_start_base,
                          new_start_base,
                          start),
           end = ifelse(end > new_end_base,
                        new_end_base,
                        end)) %>%
    mutate(width = end - start)  %>%
    dplyr::select(-new_start_base, -new_end_base)
  

  ### TO-DO: IMPLEMENT OPTION TO FILTER OUT PROBLEMATIC TRANSCRIPTS
  ### BUT KEEP AT LEAST ONE TRANSCRIPT FROM GENES THAT SEEMS TO
  ### HAVE GOOD EVIDENCE OF EXPRESSION.
    ### I essentially had this implemented in a previous version, so
    ### would just need to reimplement.
  
  ### Stringtie annotation will never have CDS, and I changed AnnotationCleaner
  ### to only operate on stringtie annotation
  # # Filter out bogus CDSs and prune ends of real CDSs
  # # Early on, I caught an edge case where an annotated CDS is partially within
  # # a tossed out exonic bin. Had to add check to filter out any CDS with a 
  # # start position downstream of new end or upstream of new start of transcript.
  # # NOTE: same filter is not needed for exons, which cannot be subject to same
  # # edge case due to which exonic bins are allowed to be filter. Can't filter 
  # # out the first exonic bin followed by a splice junction (sj) or the last exonic
  # # bin preceded by an sj. Thus, entire exons will never get filtered out of a 
  # # retained transcript.
  # CDS_annotation <- gtf %>%
  #   filter(type == "CDS") %>%
  #   inner_join(check %>%
  #                mutate(transcript_id = all_XF) %>%
  #                dplyr::select(transcript_id, new_start_base, new_end_base),
  #              by = "transcript_id") %>%
  #   filter(!(start > new_end_base) & !(end < new_start_base)) %>%
  #   group_by(transcript_id) %>%
  #   mutate(start = new_start_base,
  #          end = new_end_base) %>%
  #   mutate(width = end - start)  %>%
  #   dplyr::select(-new_start_base, -new_end_base)
  
  
  # Combine transcript, exon, and CDS information
  new_annotation <- bind_rows(list(transcript_annotation,
                                   exon_annotation)) %>%
    arrange(transcript_id)
  
  # Convert to GenomicRanges object and export
  new_GR <- GRanges(seqnames = Rle(new_annotation$seqnames),
                    ranges = IRanges(new_annotation$start, end = new_annotation$end, 
                                     names = 1:nrow(new_annotation)),
                    strand = Rle(new_annotation$strand))
  
  mcols(new_GR) <- new_annotation %>%
    dplyr::select(-seqnames, -start, -end, -strand)
  
  rtracklayer::export(new_GR,
                      con = output_path)
  
  
  
  
}



# Prune annotation -------------------------------------------------------------

### Extract Annotations used

# Original
gtf <- as_tibble(rtracklayer::import(opt$reference))

# Flattened + exon IDs
flat_gtf <- as_tibble(rtracklayer::import(opt$flatref))


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


samps <- gsub("_exonbin.featureCounts", "", files)
samps <- samps[!grepl(".summary$", samps)]
samps <- unique(samps)

# Create necessary combined table of all relevant quantifications
# for all samples.
for(i in seq_along(samps)){
  
  ### Exon bin quantification
  
  # File to extract
  filename <- paste0(opt$directory, "/", samps[i], "_exonbin.featureCounts")
  
  # Import file and modify column names
  exonbins_temp <- fread(filename)
  colnames(exonbins_temp) <- c("all_EF", "rname", "start", "end", "strand", "length", "GF", "reads")
  
  # Filter out last couple meta information rows
  exonbins_temp <- exonbins_temp %>%
    filter(!grepl("__", all_EF))
  
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

  
  ### Intron bin quantification
  
  # File to extract
  filename <- paste0(opt$directory, "/", samps[i], "_intronbin.featureCounts")
  
  # Import file and modify column names
  intronbins_temp <- fread(filename)
  colnames(intronbins_temp) <- c("all_IF", "rname", "start", "end", "strand", "length", "GF", "reads")
  
  # Filter out last couple meta information rows
  intronbins_temp <- intronbins_temp %>%
    filter(!grepl("__", all_IF))
  

  # Add mutation rate and sample ID columns
  intronbins_temp <- intronbins_temp %>%
    mutate(mutrate = 0,
           sample = samps[i]) 
  # Add it to growing final table
  if(i == 1){
    
    intronic_background  <- intronbins_temp
    
  }else{
    
    intronic_background  <- bind_rows(intronic_background, intronbins_temp)
    
  }

  
}

exonbins <- exonbins %>%
  dplyr::select(sample, GF, all_EF, reads, mutrate)

intronic_background <- intronic_background %>%
  dplyr::mutate(intron_length = length) %>%
  dplyr::select(sample, GF, all_IF, reads, mutrate, intron_length)


### Make dataframes to pass to pruner

# Map each sample to an ID corresponding to treatment condition
if(opt$bins == ""){
  
  
  ### NEED TO CHANGE THIS!!!
  if(opt$groups == ""){
    exp_id = tibble(Exp_ID = 1,
                    sample = unique(intronic_background$sample))

  }else{

    ### TO-DO: Make an option to require all exons in a group to pass
    ### Currently just need any to pass, which makes it equivalent to
    ### not having any groups
    exp_id <- as_tibble(fread(opt$groups))

  }

}else{
  
  exp_id = tibble(Exp_ID = 1,
                  sample = unique(intronic_background$sample))
  
  
}

### Score exons
dir <- dirname(opt$output)

# sample ID
if(length(samps) == 1){
  sampleID <- samps
}else{
  sampleID <- ""
}

score_exons(exonbins, flat_gtf = flat_gtf,
            gtf = gtf, dir = dir,
            exp_ids = exp_id, intronic_background = intronic_background,
            FDR = opt$fdr, coverage_floor = opt$floor, ifactor = opt$ifactor, 
            read_length = opt$readlength, tau2 = opt$priorvar,
            a1 = opt$slope, a0 = opt$intercept,
            sampleID = sampleID)

### Prune annotation

EF_to_TF_file <- paste0(dir, "/", paste0(sampleID, "EF_to_TF.csv"))
crap_exons_file <- paste0(dir, "/", paste0(sampleID, "unsupported_exons.csv"))
exon_check_file <- paste0(dir, "/", paste0(sampleID, "exon_check.csv"))

EF_to_TF <- fread(EF_to_TF_file)
crap_exons <- fread(crap_exons_file)
flat_annotation <- flat_gtf
gtf <- gtf
output_path <- opt$output
exon_check <- fread(exon_check_file)


# Create new annotation
clean_annotation(EF_to_TF = EF_to_TF,
                 crap_exons = crap_exons,
                 flat_annotation = flat_annotation,
                 gtf = gtf,
                 exon_check = exon_check,
                 dir = dir,
                 output_path = output_path)
