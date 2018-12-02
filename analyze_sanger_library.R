### Analyze sanger library

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please do not copy and/or distribute this script without proper citation of 
# multiEditR publication
###########################################################################################

### Libraries
library(sangerseqR) # To work with sanger sequencing .ab1 files
library(tidyverse) # For data manipulation
library(magrittr) # For piping
library(plyr) # For data manipulation
library(readr) # For loading and writing .txt files
source("/Users/kluesner/Desktop/Research/EditR/multiEditR/program/working_branch/dependencies.R")

##### Load files
directory = "/Users/kluesner/Desktop/Research/EditR/multiEditR/sanger_library"
files = paste0(directory, "/", list.files(directory))
ACGT = c("A", "C", "G", "T")

##### Functions
make_trimmed_df = function(file, phred_cutoff = 0.001){
  init_time = Sys.time()
  # Load sequencing files
  
  # Make sangerseq objects
  sanger = readsangerseq(file)
  
  # Generate ctrl sanger data frame
  # Generate ctrl primary basecalls
  df = make_ctrl_sanger_df(sanger)
  init_seq = df$base_call %>% paste0(., collapse = "")
  
  # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
  fastq = abif_to_fastq(path = file, cutoff = phred_cutoff)
  
  # Align the both the ctrl and samp to their fastq filtered sequences
  # reasonable to assume phred scored sequence will always be smaller than the primary seq
  # use high gap penalty to force SNP alignment
  alignment = pairwiseAlignment(pattern = fastq$seq, subject = init_seq)
  
  # Save unfiltered dataframes
  raw_df = df
  
  # Filter dfs on high phred sequence
  df %<>% 
    filter(index >= alignment@subject@range@start) %<>%
    filter(index <= alignment@subject@range@start + alignment@subject@range@width - 1) %<>%
    mutate(post_filter_index = 1:NROW(index)) %<>%
    mutate(height = {ifelse(base_call == "A", A_area,
                            ifelse(base_call == "C", C_area,
                                   ifelse(base_call == "G", G_area,
                                          ifelse(base_call == "T", T_area,NA))))}) %<>%
    filter(!grepl("N", pre_pentanucleotide)) %<>% # Filter out ends that have pre_ and post_nucleotides with NAs
    filter(!grepl("N", post_pentanucleotide)) %<>%
    mutate(path = file)
  
  # Create a smoothed regression for local peak height
  loess_model = loess(df, formula = height ~ index) 
  
  # Normalize the height based on the loess model -- heights become residuals
  df = df %>% mutate(normalized_height = height - predict(loess_model))
  
  # Convert df from wide to long about each base channel and index
  gathered_df = df %>%
    dplyr::select(everything(), -A_perc, -C_perc, -G_perc, -T_perc) %>%
    dplyr::rename(A = A_area, C = C_area, `T` = T_area, G = G_area) %>%
    gather(base, height, A:`T`) %>%
    filter(base_call == base)
  
  # Communicate the progress for each file analyzed
  stop_time = round(as.numeric(Sys.time() - init_time), 2)
  message(paste0("file ", file, " analyzed.", "\t\t\t\t", stop_time))
  
  # Return the dataframe
  return(gathered_df)
}


##### Compile data
# Takes 12.15 minutes to run through 488 files on single core
# 3.1 GHz Intel Core i5
raw_data = lapply(FUN = make_trimmed_df, X = files) %>% plyr::ldply(., data.frame)

# Write data
write_tsv(raw_data, "/Users/kluesner/Desktop/Research/EditR/multiEditR/sanger_library.tsv")

# Read data instead of writing
#raw_data = read_tsv("/Users/kluesner/Desktop/Research/EditR/multiEditR/sanger_library.tsv")

# Extract the pre nucleotide motifs
data = raw_data %>%
  mutate(pre_tetranucleotide = substr(pre_pentanucleotide, 2, 5),
         pre_trinucleotide = substr(pre_pentanucleotide, 3, 5),
         pre_dinucleotide = substr(pre_pentanucleotide, 4, 5),
         pre_nucleotide = substr(pre_pentanucleotide, 5, 5))

# Generate a pre_trinucleotide table for each motif and base of interest with median, sd and n information
prenucleotide_table = inner_join(aggregate(normalized_height ~ pre_trinucleotide*base, data, median) %>%
                                   dplyr::rename(mean = normalized_height),
                                 aggregate(normalized_height ~ pre_trinucleotide*base, data, sd) %>%
                                   dplyr::rename(sd = normalized_height) %>%
                                   inner_join(.,
                                              aggregate(normalized_height ~ pre_trinucleotide*base, data, length) %>%
                                                dplyr::rename(n = normalized_height)))

# Join the pre_trinucleotide table 
adj_data = data %>% inner_join(prenucleotide_table) %>% mutate(adj_height = normalized_height - mean)

### Plot data from all samples
alpha = 0.0001
quantiles = quantile(data$normalized_height, c((alpha/2), 1-(alpha/2)))

# Line plot
data %>%
  filter(normalized_height > quantiles[1] & normalized_height < quantiles[2]) %>%
  ggplot(., aes(x = index, y = normalized_height, color = base)) +
  #ggtitle("Unormalized Traces") +
  ylab("Normalized Height of Peak (FIU)") +
  xlab("Position of Peak") +
  labs(color = "Base") +
  scale_color_manual(values = c("A" = "forestgreen", "C" = "royalblue", "G" = "grey30", "T" = "firebrick")) +
  geom_line(alpha = 0.1) +
  geom_smooth(color = "black") +
  theme_bw(base_size = 18) +
  geom_hline(yintercept = quantiles, linetype = "dashed")


# Violin plot
data %>%
  filter(normalized_height > quantiles[1] & normalized_height < quantiles[2]) %>%
  ggplot(., aes(x = reorder(pre_trinucleotide, normalized_height), y = normalized_height, fill = base)) +
  ggtitle("Unormalized Traces") +
  ylab("Height of Peak (FIU)") +
  xlab("Position of Peak") +
  labs(color = "Base") +
  #scale_x_discrete(limits = c(-1000, 1500)) +
  scale_fill_manual(values = c("A" = "forestgreen", "C" = "royalblue", "G" = "grey30", "T" = "firebrick")) +
  geom_bar(stat = "summary", fun.y = "length", alpha = 0.5) +
  geom_violin() +
  theme_bw() +
  facet_wrap(.~base) +
  coord_flip() +
  theme(aspect.ratio = 1.2)