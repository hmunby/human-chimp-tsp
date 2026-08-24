# reformat_relate_data.R Reformat Relate Rdata files (from Speidel paper) to tab delim text file and ages file to annotate VCF with 

# age_estimate is the midpoint of the lower/upper age bounds (the point age used downstream).

args <- commandArgs(trailingOnly = TRUE)
rdata_path <- args[1]
#out_path <- args[2]
out_path_annotation <- args[2]

# Load data and export the annotation table
load(rdata_path) # loads Rdata into dataframe named allele_ages
allele_ages$age_estimate <- round((allele_ages$lower_age + allele_ages$upper_age) / 2, 1) # midpoint of the age bounds
#write.table(allele_ages, file=out_path, col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA", sep="\t")

# Reformat data to format needed to annotate VCF with allele ages 
age_annotations <- allele_ages[,c("CHR","BP","BP","lower_age","upper_age","age_estimate")]
write.table(age_annotations, file=out_path_annotation, col.names=FALSE, row.names=FALSE, quote=FALSE, na="NA", sep="\t")

q()
