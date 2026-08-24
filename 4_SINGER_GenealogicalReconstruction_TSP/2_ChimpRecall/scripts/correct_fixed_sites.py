#!/bin/usr python3 

import pandas as pd 
# ### DESCRIPTION ###
# Script to correct coding of fixed sites in .haps output from SHAPEIT

# At sites where the sample is fixed for the ALT allele during phasing, SHAPEIT recodes the REF as ALT allele and sets all GTs to 0 instead of 1

# This script corrects the coding of fixed sites in the .haps file

# For each line in file, if REF and ALT allele are different (col 4 and 5), print line 
# If REF and ALT are the same, lookup the position in the .bim file and get the REF allele, fix ref and ALT alleles and print line


haps_file = snakemake.input.haps
id_file = snakemake.input.ids
output_file = snakemake.output.haps

# # testing 
# haps_file="resources/filtered.chr1.phased.haps"
# bim_file="resources/filtered.mingeno.1.bim"
# output_file="resources/fixed.haps"

# Read in the id file to pandas df 
ids = pd.read_csv(id_file, sep='\t', header=None)
ids.columns = ['CHROM','POS','REF','ALT']

# Read haps file into pd df
haps = pd.read_csv(haps_file, sep=' ', header=None)

# For each line in haps file, check 1) If REF and ALT are same, import REF allele from ids 2) If REF and ALT are correct order as in ids 
# For discrepancies, flip 0s and 1s in cols 5 onward
for i in range(len(haps)):
    # Get the position of the site in the haps file
    chrom = haps.iloc[i,0]
    pos = haps.iloc[i,2]
    # Get the REF and ALT alleles for the site
    ref = haps.iloc[i,3]
    alt = haps.iloc[i,4]
    # Get the REF allele for the site from the ids file
    ref_id = ids[ids['POS'] == pos]['REF'].values[0]
    # If REF and ALT alleles are the same, correct the coding of the site
    if ref != ref_id:
        # Set 0s to 1s and 1s to 0s
        haps.iloc[i,5:] = 1 - haps.iloc[i,5:]
        # If REF and ALT alleles are different, correct the coding of the site
        if ref != alt:
            # Swap REF and ALT alleles
            haps.iloc[i,3] = alt
            haps.iloc[i,4] = ref
        else:
            # Set ref as ref_id
            haps.iloc[i,3] = ref_id

# Write corrected haps file to output, tab delimited
haps.to_csv(output_file, sep='\t', header=False, index=False)