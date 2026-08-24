#!/usr/bin/env Rscript

library(tidyverse)

# Import command line args
args = commandArgs(trailingOnly=TRUE)
eigenvec_file = args[1]
out_file = args[2]
sample_subspecies_file = args[3]

# Assign samples to subspecies based on eigenvals

# Eigen val cutoffs
# EAS - Pan troglodytes schweinfurthii - Eastern Chimpanzee
# WES - Pan troglodytes verus - Western Chimpanzee
# CEN - Pan troglodytes troglodytes - Central Chimpanzee

# Define cutoffs (based on known labels)
# EAS : PC1 postive & PC2 negative 
# WES : PC1 negative
# CEN : PC1 positive & PC2 positive

# Read in eigenvec file
pca <- read.table(eigenvec_file, header = FALSE)
# Remove extra ID column
pca <- pca[,-1]
# Set names
names(pca)[1] <- "Sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Add new column for subspecies based on PC1 and PC2
pca <- pca %>% 
  mutate(subspecies = case_when(
    PC1 > 0 ~ "WES",
    PC1 < 0 & PC2 < 0 ~ "CEN",
    PC1 < 0 & PC2 > 0 ~ "EAS"
  ))

# Write to file
write.table(pca, out_file, row.names = FALSE, quote = FALSE, sep = "\t")

sample_subspecies <- pca %>%
  select(Sample, subspecies)

write.table(sample_subspecies, sample_subspecies_file, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
