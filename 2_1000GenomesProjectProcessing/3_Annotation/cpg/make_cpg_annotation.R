# make_cpg_annotation.R  Generate CpG annotations for all SNPs in a VCF 

# R version 4.5.1 

# Bioconducter installations required - commented as only required upon first use
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.21")
#BiocManager::install("MutationalPatterns")
#BiocManager::install("BSgenome")
#BiocManager::install("VariantAnnotation") # Use to reload VCF and get additional info fields (ancestral states)

# Clear environment
rm(list=ls())

# Read in command-line args
args <- commandArgs(trailingOnly = TRUE)
vcf_path <- args[1]
table_out <- args[2]
# Test example: vcf_path <- "tmp/test.vcf.gz"
# Set ref genome (in BSgenome terms) - can be added as arg in future 
ref_genome <- "BSgenome.Hsapiens.1000genomes.hs37d5"

# Load libraries & reference genome 
library(MutationalPatterns)
library(VariantAnnotation)
library(stringr)
library(ref_genome, character.only = TRUE)

# Read VCF -- must already be only SNPs 
vcf <- read_vcfs_as_granges(vcf_path, NA, ref_genome, predefined_dbs_mbs = TRUE)
vcf <- vcf[[1]]

# Read info 
info_AA <- readInfo(vcf_path, "AA", param=ScanVcfParam())
# Keep only rows where name is in vcf object (VCF has duplicate positions/multiallelic sites removed)
info_AA <- info_AA[which(names(info_AA) %in% names(vcf))]
# Remove duplicates (keep first instance)
info_AA <- info_AA[!duplicated(names(info_AA))]
# Add to vcf object
vcf$AA <- info_AA
# Reformat AA to ancestral allele only and clean up unknowns 
vcf$AA <- str_split_fixed(vcf$AA, "", 2)[,1]
vcf$AA[vcf$AA == "."] <- "N"
vcf$AA[vcf$AA == "-"] <- "N"
# Capitalise
vcf$AA <- toupper(vcf$AA)

# Mut type context 
type_context <- type_context(vcf, ref_genome)
type_context <- data.frame(type_context$types, as.character(type_context$context))
colnames(type_context) <- c("type","context")
# For testing purposes add other columns to type_context (not required for applying the function but here for clarity)
type_context$AA <- vcf$AA
# Ref and alt , convert to character first 
type_context$REF <- as.character(vcf$REF)
type_context$ALT <- as.character(unlist(vcf$ALT))

# Define function to determine whether a mutation is a CpG based on mutation type and 3n context based on REF alone (not used for final annotation, but was used to calculate error rates when comparing to ancestral allele-based annotation)
is.CpG <- function(mut, con){
  con_cg = str_sub(con,2,-1)
  if(mut == "C>T" & con_cg == "CG"){
    cpg = 1
  }else{
    cpg = 0
  }
  return(cpg)
}

# Define function to determine whether a mutation is a CpG based on mutation type and 3n context based on AA. Ancestral allele info is used to identify CpGs in cases where ref/alt do not reflect ancestral state. First statement takes care of false positive errors (i.e. it is actually a T->C in a TG context) and second statement takes care of false negative errors (i.e. it is actually a C->T in a CG context)
is.CpG.AA <- function(mut, con, AA, REF, ALT){
  con_cg = str_sub(con,2,-1)
  if(AA == "N"){
    # Use original criterion if ancestral allele is unknown
    if(mut == "C>T" & con_cg == "CG"){
      cpg = 1
    }else{
      cpg = 0
    }
  }else{
    if(mut == "C>T" & con_cg == "CG" & AA == REF){
      cpg = 1
    }else if(mut == "T>C" & con_cg == "TG" & AA == ALT){
      cpg = 1
    }else{
      cpg = 0
    }
  }
  return(cpg)
}

# Add CpG column to vcf object, remove type_context object
vcf$CpG <- mapply(is.CpG, type_context$type, type_context$context)
# Use extended function with ancestral allele info to refine CpG annotations
vcf$CpG_AA <- mapply(is.CpG.AA, type_context$type, type_context$context, type_context$AA, type_context$REF, type_context$ALT)

# Compare CpG and CpG_AA columns
type_1_error_rate <- sum(vcf$CpG == 1 & vcf$CpG_AA == 0) / sum(vcf$CpG == 1)
type_2_error_rate <- sum(vcf$CpG == 0 & vcf$CpG_AA == 1) / sum(vcf$CpG == 0)
print(paste("Type 1 error rate (false positive CpG annotation):", round(type_1_error_rate,4)))
print(paste("Type 2 error rate (false negative CpG annotation):", round(type_2_error_rate,4)))

# Final % CpG 
print(paste("Final % CpG SNPs:", round((sum(vcf$CpG_AA == 1) / length(vcf$CpG_AA))*100,2)))

rm(type_context)
# Use decode to build data.frame() of columns to be written out 
vcf <- data.frame(decode(vcf@seqnames), (vcf@ranges@start), vcf$REF, unlist(vcf$ALT), vcf$CpG_AA)
write.table(vcf, file = table_out ,row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")

quit()