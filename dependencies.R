### Dependencies for multiEditR
### Measure GC Content by exon of Ensembl transcript IDs

#############################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
#############################################################################

### Libraries
library(sangerseqR)
library(tidyverse)
library(magrittr)
library(plyr)
library(sangerseqR)
library(gamlss)
library(CrispRVariants)
library(readr)

### Data
bases = c("A", "C", "G", "T")

pre_trinucleotide_table = read_tsv("/Users/kluesner/Desktop/Research/EditR/multiEditR/pre_trinucleotide_table.tsv")

phred_scores = data.frame(stringsAsFactors=FALSE,
                          phred = c("!", "“", "$", "%", "&", "‘", "(", ")", "*", "+", ",", "–",
                                    ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                    ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F",
                                    "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                                    "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_",
                                    "`", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
                                    "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y",
                                    "z", "{", "|", "}", "~"),
                          prob = c(1, 0.794328235, 0.501187234, 0.398107171, 0.316227766,
                                   0.251188643, 0.199526232, 0.158489319, 0.125892541, 0.1,
                                   0.079432824, 0.063095734, 0.050118723, 0.039810717, 0.031622777,
                                   0.025118864, 0.019952623, 0.015848932, 0.012589254, 0.01,
                                   0.007943282, 0.006309573, 0.005011872, 0.003981072, 0.003162278,
                                   0.002511886, 0.001995262, 0.001584893, 0.001258925, 0.001, 0.000794328,
                                   0.000630957, 0.000501187, 0.000398107, 0.000316228, 0.000251189,
                                   0.000199526, 0.000158489, 0.000125893, 1e-04, 7.94328e-05,
                                   6.30957e-05, 5.01187e-05, 3.98107e-05, 3.16228e-05, 2.51189e-05,
                                   1.99526e-05, 1.58489e-05, 1.25893e-05, 1e-05, 7.9433e-06,
                                   6.3096e-06, 5.0119e-06, 3.9811e-06, 3.1623e-06, 2.5119e-06, 1.9953e-06,
                                   1.5849e-06, 1.2589e-06, 1e-06, 7.943e-07, 6.31e-07, 5.012e-07,
                                   3.981e-07, 3.162e-07, 2.512e-07, 1.995e-07, 1.585e-07, 1.259e-07,
                                   1e-07, 7.94e-08, 6.31e-08, 5.01e-08, 3.98e-08, 3.16e-08,
                                   2.51e-08, 2e-08, 1.58e-08, 1.26e-08, 1e-08, 7.9e-09, 6.3e-09, 5e-09,
                                   4e-09, 3.2e-09, 2.5e-09, 2e-09, 1.6e-09, 1.3e-09, 1e-09, 8e-10,
                                   6e-10, 5e-10)
)

### Functions
# Function to analyze an individual .ab1 file
# Need to make sure the file ends in .ab1 for this function
make_ctrl_sanger_df = function(sanger_file){
  base_calls = makeBaseCalls(sanger_file)
  sanger_df = base_calls %>% peakAmpMatrix %>% data.frame()
  colnames(sanger_df) = c("A_area","C_area","G_area","T_area")
  sanger_df %<>% 
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = 100*A_area / Tot.Area,
           C_perc = 100*C_area / Tot.Area,
           G_perc = 100*G_area / Tot.Area,
           T_perc = 100*T_area / Tot.Area) %<>%
    mutate(base_call = strsplit(x = toString(base_calls@primarySeq), split = "") %>% unlist) %<>%
    mutate(index = 1:NROW(.)) %<>%
    mutate(pre_trinucleotide = paste0(lag(base_call, n = 3),
                                      lag(base_call, n = 2),
                                      lag(base_call, n = 1))) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                            ifelse(max_base == "C", C_area,
                                   ifelse(max_base == "G", G_area,
                                          ifelse(max_base == "T", T_area,NA))))}) 
}


### Makes a df for the edited sample that accounts for the misalignment of peaks
### Employs the secondary basecalls
### Not aligned for peaks that are shifted by a complete basecall
make_samp_sanger_df = function(samp_sanger, ctrl_seq){
  
  ### Align the phase of the primary and secondary basecalls in the sample sequence to that of the control
  ### This changes where and what the bases are called as
  ### Returns an object of class sangerseq
  phased_samp_sanger = samp_sanger %>%
    makeBaseCalls(.) %>%
    setAllelePhase(., ctrl_seq)
  
  ### Return a data frame from the phased sample sangerseq object with the position of each 
  ### This method appears to return higher intensities for the noise, thus we're going to trust it more for detecting noise for modelling
  samp_peakAmpDF = phased_samp_sanger@peakPosMatrix %>%
    as.data.frame() %>%
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    na_if(., 0) %>%
    dplyr::mutate(primary_base_call = primarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist(),
                  secondary_base_call = secondarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist()) %>%
    dplyr::mutate(identical = {primary_base_call == secondary_base_call}) %>%
    dplyr::mutate(row_sum = as.integer(!is.na(A)) +
                    as.integer(!is.na(C)) +
                    as.integer(!is.na(G)) +
                    as.integer(!is.na(`T`))
    ) %>% 
    dplyr::mutate(peak_pos = {ifelse(grepl("A|C|G|T", primary_base_call),
                                     ifelse(primary_base_call == "A", A, 
                                            ifelse(primary_base_call == "C", C, 
                                                   ifelse(primary_base_call == "G", G, 
                                                          ifelse(primary_base_call == "T", `T`, "Error, stop!")))),
                                     pmin(A, C, G, `T`, na.rm = TRUE))}) %>%
    dplyr::mutate(A = ifelse(is.na(A), peak_pos, A) %>% as.numeric(),
                  C = ifelse(is.na(C), peak_pos, C) %>% as.numeric(),
                  G = ifelse(is.na(G), peak_pos, G) %>% as.numeric(),
                  `T` = ifelse(is.na(`T`), peak_pos, `T`) %>% as.numeric())
  
  ### Reformat df to be identical to output of make_ctrl_sanger_df()
  samp_peakAmpDF %<>%
    dplyr::mutate(A_area = samp_sanger@traceMatrix[samp_peakAmpDF$A, 1],
                  C_area = samp_sanger@traceMatrix[samp_peakAmpDF$C, 2],
                  G_area = samp_sanger@traceMatrix[samp_peakAmpDF$G, 3],
                  T_area = samp_sanger@traceMatrix[samp_peakAmpDF$`T`, 4]) %<>%
    dplyr::select(A_area:T_area, primary_base_call, secondary_base_call) %<>%
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = 100*A_area / Tot.Area,
           C_perc = 100*C_area / Tot.Area,
           G_perc = 100*G_area / Tot.Area,
           T_perc = 100*T_area / Tot.Area) %>%
    mutate(index = 1:NROW(Tot.Area)) %<>%
    mutate(pre_trinucleotide = paste0(lag(primary_base_call, n = 3),
                                      lag(primary_base_call, n = 2),
                                      lag(primary_base_call, n = 1))) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                                     ifelse(max_base == "C", C_area,
                                            ifelse(max_base == "G", G_area,
                                                   ifelse(max_base == "T", T_area,NA))))}) 
  
  return(samp_peakAmpDF)
}

# The alignment index and subsequent filtering is currently off -- need to adapt to take the longest region of alignment, and then only filter our those sequences.
align_sanger_dfs = function(control_df, sample_df){
  
  ctrl_seq = control_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  sample_seq = sample_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  
  if(ctrl_seq@length > sample_seq@length){
    alignment = pairwiseAlignment(pattern = sample_seq, subject = ctrl_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
  } else
  {
    alignment = pairwiseAlignment(pattern = ctrl_seq, subject = sample_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
  }
  return(list(control_df, sample_df))
}

### Convert the abif to fastq
abif_to_fastq = function (seqname = "sample", path, trim = TRUE, cutoff = 1, 
                          min_seq_len = 20, offset = 33, recall = FALSE) 
{
  sangerseqr <- requireNamespace("sangerseqR")
  stopifnot(isTRUE(sangerseqr))
  abif <- sangerseqR::read.abif(path)
  if (is.null(abif@data$PCON.2)) {
    message(sprintf("failed on %s", seqname))
    return()
  }
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
  if (!typeof(abif@data$PCON.2) == "integer") {
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)]
  }
  else {
    num_quals <- abif@data$PCON.2[1:length(abif@data$PLOC.2)]
  }
  if (isTRUE(recall)) {
    recalled <- sangerseqR::makeBaseCalls(sangerseqR::sangerseq(abif))
    nucseq <- sangerseqR::primarySeq(recalled, string = TRUE)
    if (nchar(nucseq) != length(num_quals)) {
      trim <- FALSE
      num_quals <- rep(60, nchar(nucseq))
      warning("Length of quality scores does not equal length of\n              re-called base sequence, ignoring quality scores")
    }
  }
  if (trim == FALSE) {
    tmp1 = list(seqname = seqname, seq = nucseq, 
                quals = rawToChar(as.raw(num_quals + offset)))
    return(tmp1)
  }
  trim_msg <- "Sequence %s can not be trimmed because it is shorter than the trim\n               segment size"
  if (nchar(nucseq) <= min_seq_len) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  scores = cutoff - 10^(num_quals/-10)
  running_sum <- rep(0, length(scores) + 1)
  for (i in 1:length(scores)) {
    num <- scores[i] + running_sum[i]
    running_sum[i + 1] <- ifelse(num < 0, 0, num)
  }
  trim_start <- min(which(running_sum > 0)) - 1
  trim_finish <- which.max(running_sum) - 2
  if (trim_finish - trim_start < min_seq_len - 1) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  tmp2 = list(seqname = seqname, seq = substring(nucseq, 
                                                 trim_start, trim_finish), quals = rawToChar(as.raw(num_quals[trim_start:trim_finish] + 
                                                                                                      offset)))
  return(tmp2)
}



trim_ends = function(fastq, threshold){
  
}

### Takes sequence and runs through a while loop until they are properly trimmed
### Appears for CTSS 11 sample
align_and_trim = function(pattern_seq, subject_seq, min_continuity = 15){
  
  gap_length = min_continuity - 1
  raw_pattern = lapply(FUN = rep, X = 1:gap_length, x = "[A-Z]") %>%
    lapply(FUN = paste0, X = ., collapse = "") %>%
    unlist()
  gsub_pattern = raw_pattern %>%
    paste0("^", ., "-|") %>%
    paste0(collapse = "") %>%
    paste0(., raw_pattern %>%
             paste0("-", ., "$|") %>%
             paste0(collapse = ""),
           collapse = "") %>%
    gsub('.{1}$','', .)
  
  input_pattern = pattern_seq
  input_subject = subject_seq
  
  alignment = pairwiseAlignment(pattern = pattern_seq, subject = subject_seq)
  
  output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  
  while(input_pattern != output_pattern | input_subject != output_subject){
    alignment = pairwiseAlignment(pattern = output_pattern, subject = output_subject)
    
    old_output_pattern = output_pattern
    old_output_subject = output_subject
    
    new_output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    new_output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    
    output_pattern = new_output_pattern
    output_subject = new_output_subject
    
    input_pattern = old_output_pattern
    input_subject = old_output_subject
  }
  
  return(list("pattern" = alignment@pattern %>% as.character() %>% gsub("-", "", .),
              "subject" = alignment@subject %>% as.character() %>% gsub("-", "", .),
              "alignment" = alignment))
}

### Substitute multiple string positions simultaneously
subchar <- function(string, pos, char) { 
  for(i in pos) { 
    substr(string, i, i) = char
  } 
  string 
} 

