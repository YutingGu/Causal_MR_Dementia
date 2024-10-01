#
# example function for IV selection
# - read summary-level genetic associations
# - filter by p-value 
# - align with reference file and  
# - clump (remove correlated genetic variants)
#
setwd("/rds/general/project/hda_23-24/live/TDS/group1/Scripts")
# packages
library("plyr") 
library("tidyr") 
#library("tidyverse")
library("stringr") 
library("vroom") 
library("data.table")  
library("ieugwasr")

#
# https://www.ebi.ac.uk/gwas/studies/GCST90274812
# PUBMED https://www.ebi.ac.uk/gwas/publications/37563310
#

# read in the data and rename column names
IL33.data = vroom("../GWAS_Data/GCST90274812.tsv.gz")
ST2.data = vroom("../GWAS_Data/ST2.txt.gz")


#chr (4): effect_allele, other_allele, variant_id, rsid
#dbl (7): chromosome, base_pair_location, beta, standard_error, effect_allele...
IL33.data = IL33.data %>% dplyr::select("rsid","effect_allele","other_allele","beta","standard_error","p_value")
IL33.data = IL33.data %>% dplyr::rename(beta_not_aligned = beta, se=standard_error, pval = "p_value", EA=effect_allele, NEA=other_allele)
dim(IL33.data)
n_distinct(IL33.data$rsid) # number of distinct SNPs

ST2.data = ST2.data %>% dplyr::select("MarkerName","Allele1","Allele2","Effect","StdErr","P-value")
ST2.data = ST2.data %>% dplyr::rename(beta_not_aligned = Effect, se=StdErr, pval = "P-value", EA=Allele1, NEA=Allele2)
ST2.data$MarkerName = gsub("_",":", ST2.data$MarkerName) # organise string format
ST2.data$EA = toupper(ST2.data$EA) # organise string format
ST2.data$NEA = toupper(ST2.data$NEA) # organise string format
ST2.data = ST2.data %>% dplyr::rename(variant = MarkerName)
dim(ST2.data)
n_distinct(ST2.data$variant)
#
# Read in the annotation file README https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291
# https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz
#
annot = vroom("../GWAS_Data/variants.tsv.bgz")
annot = annot %>% dplyr::select("variant","chr","pos","rsid", "ref", "alt")
dim(annot)
n_distinct(annot$rsid)
n_distinct(annot$variant)
#annot = annot %>% unite("MarkerName", chr:pos, sep = ":")

#pos	Position of the variant in GRCh37 coordinates.
#ref    string  Reference allele on the forward strand.
#alt    string  Alternate allele (not necessarily minor allele)
#beta   float   Estimated effect size of alt allele.

IL33.data = merge(x=IL33.data,y=annot,by="rsid",sort=FALSE)
dim(IL33.data)
n_distinct(IL33.data$rsid)
head(IL33.data)


ST2.data = merge(x=ST2.data,y=annot,by="variant",sort=FALSE)
dim(ST2.data)
n_distinct(ST2.data$variant)
rm(annot)

#
# Align the effect allele according to the annotation file
# IL33
table(IL33.data$alt == IL33.data$EA)
table(IL33.data$alt == IL33.data$NEA)
#odd matches
head(IL33.data[which(IL33.data$alt != IL33.data$EA & IL33.data$alt != IL33.data$NEA),])
#remove those
IL33.data = IL33.data[-which(IL33.data$alt != IL33.data$EA & IL33.data$alt != IL33.data$NEA),]
dim(IL33.data)
n_distinct(IL33.data$rsid)

table(IL33.data$alt == IL33.data$EA)
table(IL33.data$alt == IL33.data$NEA)
# align beta value
IL33.data$beta = ifelse(IL33.data$alt==IL33.data$EA, IL33.data$beta_not_aligned, -1*IL33.data$beta_not_aligned)

################
# ST2
# run code below to check there is no odd matches
table(ST2.data$alt == ST2.data$EA)
table(ST2.data$alt == ST2.data$NEA)
# align beta value
ST2.data$beta = ifelse(ST2.data$alt==ST2.data$EA, ST2.data$beta_not_aligned, -1*ST2.data$beta_not_aligned)


#
# Filter by p-value (please note there are no genetic variants below genome-wide significance, I have lowered the p-value)
#
IL33_iv_out = IL33.data[IL33.data$pval<10^-6,]
dim(IL33_iv_out)
ST2_iv_out = ST2.data[ST2.data$pval<10^-6,]
dim(ST2_iv_out)
ST2_iv_out2 = ST2.data[ST2.data$pval<5e10^-8,]
dim(ST2_iv_out2)
#
# Clumping at different r2 thresholds (r2 - 0.01 or 0.001 is recommended)
#

#data_clump = ieugwasr::ld_clump(iv_out)
#dim(data_clump)
IL33_clump_01 = ieugwasr::ld_clump(IL33_iv_out, clump_r2 = 0.01)
dim(IL33_clump_01)
ST2_clump_01 = ieugwasr::ld_clump(ST2_iv_out, clump_r2 = 0.01)
dim(ST2_clump_01)
ST2_clump_02 = ieugwasr::ld_clump(ST2_iv_out2, clump_r2 = 0.01)
dim(ST2_clump_02)

################
# load the saved data
save(IL33_iv_out, IL33_clump_01, ST2_iv_out, ST2_clump_01, ST2_clump_02, file = "../GWAS_Data/iv_out.RData") 
load(file = "../GWAS_Data/iv_out.RData")

quit(save="no")


# R CMD BATCH 1-read-filter-clump-group1.R 1.out &


