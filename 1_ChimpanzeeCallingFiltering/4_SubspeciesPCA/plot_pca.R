#!/usr/bin/env Rscript
# Supplementary Figure 1: PC1 vs PC2 of the LD-pruned chimpanzee genotypes.
#
# Points are coloured by the subspecies cluster they fall in and shaped by whether that subspecies
# was already known from the sample's source publication (circle) or assigned here (triangle).
# Axis labels carry the percentage of variance explained.
#
# Usage: plot_pca.R <pca.eigenvec> <pca.eigenval> <sample_known_subspecies.txt> <out.pdf>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: plot_pca.R <eigenvec> <eigenval> <known_subspecies> <out.pdf>")
}
eigenvec_file <- args[1]
eigenval_file <- args[2]
known_subspecies_file <- args[3]
out_file <- args[4]

# plink writes FID and IID as the first two columns; --double-id makes them identical, so drop one.
pca <- read.table(eigenvec_file, header = FALSE)
pca <- pca[, -1]
names(pca)[1] <- "Sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))
pca$Sample <- as.character(pca$Sample)

eigenval <- scan(eigenval_file)
pve <- 100 * eigenval / sum(eigenval)

# Cluster -> subspecies. The boundaries come from where the already-labelled samples fall:
#   WES (P. t. verus)          PC1 > 0
#   CEN (P. t. troglodytes)    PC1 < 0, PC2 < 0
#   EAS (P. t. schweinfurthii) PC1 < 0, PC2 > 0
pca <- pca %>%
  mutate(subspecies = case_when(
    PC1 > 0              ~ "WES",
    PC1 < 0 & PC2 < 0    ~ "CEN",
    PC1 < 0 & PC2 > 0    ~ "EAS"
  ))

# "UNK" in the known-subspecies file marks the samples with no previously published label.
known <- read.table(known_subspecies_file, header = FALSE,
                    col.names = c("Sample", "known_label"),
                    colClasses = "character")
known <- known %>%
  mutate(label_source = ifelse(known_label != "UNK", "Known", "Unknown"))

pca <- pca %>% left_join(known, by = "Sample")

colors <- c("WES" = "purple", "CEN" = "red", "EAS" = "blue")
shapes <- c("Known" = 16, "Unknown" = 17)

p <- ggplot(pca, aes(x = PC1, y = PC2, color = subspecies, shape = label_source)) +
  geom_point(size = 2.5) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]),
       y = sprintf("PC2 (%.1f%%)", pve[2]),
       color = "Subspecies", shape = NULL) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_classic()

ggsave(out_file, p, width = 8, height = 8)
