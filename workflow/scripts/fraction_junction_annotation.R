### PURPOSE OF THIS SCRIPT
## Determine what fraction of high quality junctions in your alignment
## are present in the annotation

# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(dplyr)
library(data.table)


### FILE PATHS

# Path to STAR's SJ.out.tab file
sj_out <- "G:/Shared drives/Matthew_Simon/IWV/Hogg_lab/fastq2EZbakR/jm_ensembl_stringtie/4hr/align/11j_4hr_1-SJ.out.tab"

# Path to annotation to score
gtf_path <- "G:/Shared drives/Matthew_Simon/IWV/Annotations/Hogg_annotations/Simple_Annotations/Hs_ensembl_lvl1_stringtie.gtf"

### PARAMETERS FOR FILTERING STAR JUNCTIONS
# Mimics what can be done by STAR automatically with following parameters:
# 1) --outSJfilterCountUniqueMin
# 2) --outSJfilterCountTotalMin
# Also filters for smallest allowed max overhang (kinda like --outSJfilterOverhangMin)
# Defaults essentially don't filter the SJ.out.tab file any more than
# STAR already filtered it.

# Smallest max overhang allowed
  # Vector of length 4th
  # 1st element: Cutoff for non-canonical intron motifs
  # 2nd element: Cutoff for GT/AG and CT/AC motifs
  # 3rd element: Cutoff for GC/AG and CT/GC motifs
  # 4th element: Cutoff for AT/AC and GT/AT motifs
min_maxoverhang <- c(1000, 1000, 1000, 1000)

# Minimum number of uniquely mapping reads supporting junction
  # Vector of length 4th
  # 1st element: Cutoff for non-canonical intron motifs
  # 2nd element: Cutoff for GT/AG and CT/AC motifs
  # 3rd element: Cutoff for GC/AG and CT/GC motifs
  # 4th element: Cutoff for AT/AC and GT/AT motifs
min_unique <- c(3, 1, 1, 1)

# Minimum number of total reads supporting junction
  # Vector of length 4th
  # 1st element: Cutoff for non-canonical intron motifs
  # 2nd element: Cutoff for GT/AG and CT/AC motifs
  # 3rd element: Cutoff for GC/AG and CT/GC motifs
  # 4th element: Cutoff for AT/AC and GT/AT motifs
min_tot <- c(3, 1, 1, 1)


# Are junctions discovered by STAR in annotation? ------------------------------


# Load in splice junctions
sjs <- read.delim(file.path(sj_out),
                  header = FALSE, as.is = TRUE) %>%
  setNames(c("seqnames", "start", "end", "strand", "motif", "annot",
             "uniqreads", "mmreads", "maxoverhang")) %>% 
  dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>%
  dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>%
  dplyr::mutate(seqnames = as.character(seqnames))
sjs <- sjs %>%
  as_tibble() %>%
  dplyr::filter(strand != "0") %>%
  dplyr::mutate(trash = case_when(
    motif == 0 & 
      (uniqreads > min_unique[1]) & 
      ((uniqreads + mmreads) > min_tot[1]) ~ FALSE,
    motif == 1 & 
      (uniqreads > min_unique[2]) & 
      ((uniqreads + mmreads) > min_tot[2]) ~ FALSE,
    motif == 2 & 
      (uniqreads > min_unique[2]) & 
      ((uniqreads + mmreads) > min_tot[2]) ~ FALSE,
    motif == 3 & 
      (uniqreads > min_unique[3]) & 
      ((uniqreads + mmreads) > min_tot[3]) ~ FALSE,
    motif == 4 & 
      (uniqreads > min_unique[3]) & 
      ((uniqreads + mmreads) > min_tot[3]) ~ FALSE,
    motif == 5 & 
      (uniqreads > min_unique[4]) & 
      ((uniqreads + mmreads) > min_tot[4]) ~ FALSE,
    motif == 6 & 
      (uniqreads > min_unique[4]) & 
      ((uniqreads + mmreads) > min_tot[4]) ~ FALSE,
    .default = TRUE
  )
  ) %>%
  dplyr::filter(!trash) %>%
  dplyr::select(-trash)
  

# Load in annotation
gtf <- rtracklayer::import(file.path(gtf_path))
annotated_sjs <- gtf %>%
  as_tibble() %>%
  dplyr::filter(type == "exon") %>%
  dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, exon_number) %>%
  dplyr::distinct() %>%
  group_by(gene_id, transcript_id) %>%
  dplyr::arrange(gene_id, transcript_id, exon_number) %>%
  dplyr::mutate(junction_start = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = end + 1
  ),
  junction_end = case_when(
    exon_number == max(exon_number) ~ NA,
    .default = lead(start, order = exon_number) - 1
  )) %>%
  na.omit() %>%
  dplyr::ungroup() %>%
  dplyr::select(seqnames, junction_start, junction_end,
                strand) %>%
  dplyr::rename(
    start = junction_start,
    end = junction_end
  )

# How many splice junctions are in both?
sj_in_both <- sjs %>%
  inner_join(
    annotated_sjs,
    by = c("seqnames", "start", "end", "strand")
  )


# Calculate score
number_of_star_junctions <- nrow(sjs)
number_in_annotation <- nrow(sj_in_both)

fraction_junction_score <- number_in_annotation / number_of_star_junctions
fraction_junction_score
