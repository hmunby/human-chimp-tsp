#!/bin/usr python3 

import pandas as pd 
import sys
import numpy as np

# Import from stdin 
samples_file=sys.argv[1]
sites_file=sys.argv[2]
shared_pos_file=sys.argv[3]
n=int(sys.argv[4])
output_file=sys.argv[5]

# Test args 
# samples_file="test/test.samples"
# sites_file="test/test.sites"
# shared_pos_file="misc_files/top_200_mean_rank.snp_pos.txt"
# n=50
# output_file="test/test_downsampled.keep"

# Read samples file into list
samples = pd.read_csv(samples_file, header=None).values.flatten()
# Add CHROM and POS to samples list to make header
samples = ["CHROM", "POS"] + list(samples)

# Read sites file into pd df
sites = pd.read_csv(sites_file, sep='\t', header=None)
sites.columns = samples

# Read shared_pos_file into pandas df 
shared_pos = pd.read_csv(shared_pos_file, sep='\t', header=0)
shared_pos.columns = ["CHROM", "POS"]

# Subset sites to any site in shared_pos (matching CHROM and POS)
sites = sites[sites["POS"].isin(shared_pos["POS"])]

# If sites is now empty raise error 
if sites.empty:
    raise ValueError("No shared sites between sites and shared_pos files")

# Keep only first line of sites
sites = sites.head(1)

# Recode genotypes to numbers - "1|0" and "0|1" to 1, "0|0" to 0, "1|1" to 2
sites = sites.replace("1|0", 1)
sites = sites.replace("0|1", 1)
sites = sites.replace("0|0", 0)
sites = sites.replace("1|1", 2)

# Subset to human - cols starting with "HG" or "NA"
human_cols = [col for col in sites.columns if col.startswith("HG") or col.startswith("NA")]
# Len human cols 
n_human = len(human_cols)
sites_human = sites[["CHROM", "POS"] + human_cols]

# For each genotype class (0,1,2), get col IDs, and take a random sample of n/freq(class)

genotype_0 = np.array(sites_human.columns[sites_human.iloc[0] == 0])
genotype_0_n = n*len(genotype_0)/n_human

genotype_1 = np.array(sites_human.columns[sites_human.iloc[0] == 1])
genotype_1_n = n*len(genotype_1)/n_human

genotype_2 = np.array(sites_human.columns[sites_human.iloc[0] == 2])
genotype_2_n = n*len(genotype_2)/n_human

# Work out how to round up/down to nearest whole number to get n samples
if round(genotype_0_n) + round(genotype_1_n) + round(genotype_2_n) == (n+1):
    if round(genotype_0_n) > 1:
        genotype_0_n = round(genotype_0_n) - 1
    else:
        genotype_1_n = round(genotype_1_n) - 1
elif round(genotype_0_n) + round(genotype_1_n) + round(genotype_2_n) == (n-1):
    genotype_0_n = round(genotype_0_n) + 1

downsample_0 = np.random.choice(genotype_0, round(genotype_0_n), replace=False)
downsample_1 = np.random.choice(genotype_1, round(genotype_1_n), replace=False)
downsample_2 =  np.random.choice(genotype_2, round(genotype_2_n), replace=False)

# Combine downsampled cols
downsampled_cols = np.concatenate([downsample_0, downsample_1, downsample_2])

# Order downsamples cols so they are in same order as original sites df - alphabetic 
downsampled_cols = np.sort(downsampled_cols)

# List of chimp samples - not human and not CHROM or POS
chimp_samples = [col for col in sites.columns if col not in human_cols]
chimp_samples = [col for col in chimp_samples if col not in ["CHROM", "POS"]]

downsample_list = list(chimp_samples) + list(downsampled_cols) 

# Output to file, one sample ID per line, no quotes
with open(output_file, "w") as f:
    for sample in downsample_list:
        f.write(sample + "\n")


