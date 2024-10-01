library("plyr") 
library("tidyverse")
library("stringr") 
library("vroom") 
library("data.table") 
library("MendelianRandomization") 
library("ieugwasr")

#################################################################################
# 1. Read in the data for two proteins and create the genetic instrument at genome wide signifcance level 
#################################################################################
# load protein summary GWAS data
IL33.data = vroom("../GWAS_Data/GCST90274812.tsv.gz")
ST2.data = vroom("../GWAS_Data/ST2.txt.gz")

# select relevant columns and rename for convenience
IL33.data = IL33.data %>% dplyr::select("rsid","effect_allele","other_allele","beta","standard_error","p_value")
IL33.data = IL33.data %>% dplyr::rename(beta_not_aligned = beta, se=standard_error, pval = "p_value", EA=effect_allele, NEA=other_allele)

ST2.data = ST2.data %>% dplyr::select("MarkerName","Allele1","Allele2","Effect","StdErr","P-value")
ST2.data = ST2.data %>% dplyr::rename(beta_not_aligned = Effect, se=StdErr, pval = "P-value", EA=Allele1, NEA=Allele2)
ST2.data$MarkerName = gsub("_",":", ST2.data$MarkerName) # organise string format
ST2.data$EA = toupper(ST2.data$EA) # organise string format
ST2.data$NEA = toupper(ST2.data$NEA) # organise string format
ST2.data = ST2.data %>% dplyr::rename(variant = MarkerName)

# load UK Biobank annotation file
annot = vroom("../GWAS_Data/variants.tsv.bgz")
annot = annot %>% dplyr::select("variant","chr","pos","rsid", "ref", "alt")

# merge UK biobank annotition file with IL33/ST2 summary GWAS data
IL33.data = merge(x=IL33.data,y=annot,by="rsid",sort=FALSE)
ST2.data = merge(x=ST2.data,y=annot,by="variant",sort=FALSE)

# filter instrument variables above GWAS significance level
IL33.data.iv = IL33.data[which(IL33.data$pval <10^-6),]
ST2.data.iv1 = ST2.data[which(ST2.data$pval <10^-6),]
ST2.data.iv2 = ST2.data[which(ST2.data$pval <5e10^-8),]

#################################################################################
#2. Read in the data for Alzheimer's disease
#################################################################################
# load AD summary GWAS data
AD.data <- vroom("../GWAS_Data/AD_GWAS.txt")

AD.data = AD.data %>% select("MarkerName","Effect_allele","Non_Effect_allele","Beta","SE","Pvalue")
AD.data = AD.data %>% dplyr::rename(rsid = MarkerName,
                                    AD_effectallele = Effect_allele, 
                                    AD_noneffectallele = Non_Effect_allele, 
                                    AD_beta_not_aligned = Beta, 
                                    AD_se = SE, 
                                    AD_pval = Pvalue)
#################################################################################
# 3. Merge the two data sets
#################################################################################
IL33.merge_AD =merge(IL33.data.iv, AD.data, by=c("rsid"))
ST2.merge_AD1 =merge(ST2.data.iv1, AD.data, by=c("rsid"))
ST2.merge_AD2 =merge(ST2.data.iv2, AD.data, by=c("rsid"))

#################################################################################
# 4. Harmonise and Prune the merged data 
#################################################################################
# Create a TRUE/FALSE vector to indicate if EA is equal to AD_effectallele. 
dim(IL33.merge_AD)
table(toupper(IL33.merge_AD$EA) == IL33.merge_AD$AD_effectallele)
table(toupper(IL33.merge_AD$EA) == IL33.merge_AD$AD_noneffectallele)

dim(ST2.merge_AD1)
table(toupper(ST2.merge_AD1$EA) == ST2.merge_AD1$AD_effectallele)
table(toupper(ST2.merge_AD1$EA) == ST2.merge_AD1$AD_noneffectallele)

dim(ST2.merge_AD2)
table(toupper(ST2.merge_AD2$EA) == ST2.merge_AD2$AD_effectallele)
table(toupper(ST2.merge_AD2$EA) == ST2.merge_AD2$AD_noneffectallele)

# The harmonisation of the effect alleles can be performed as below in this case we are just creating the variable AD_beta as all the alleles are aligned 
IL33.merge_AD$AD_beta = ifelse((IL33.merge_AD$EA) == IL33.merge_AD$AD_effectallele, 
                               IL33.merge_AD$AD_beta_not_aligned, -1*IL33.merge_AD$AD_beta_not_aligned)
ST2.merge_AD1$AD_beta = ifelse((ST2.merge_AD1$EA) == ST2.merge_AD1$AD_effectallele, 
                               ST2.merge_AD1$AD_beta_not_aligned, -1*ST2.merge_AD1$AD_beta_not_aligned)
ST2.merge_AD2$AD_beta = ifelse((ST2.merge_AD2$EA) == ST2.merge_AD2$AD_effectallele, 
                               ST2.merge_AD2$AD_beta_not_aligned, -1*ST2.merge_AD2$AD_beta_not_aligned)

table(IL33.merge_AD$AD_beta == IL33.merge_AD$AD_beta_not_aligned)
table(ST2.merge_AD1$AD_beta == ST2.merge_AD1$AD_beta_not_aligned)
table(ST2.merge_AD2$AD_beta == ST2.merge_AD2$AD_beta_not_aligned)


#Clump and prune

#The final step is to prune or clump the SNPs. Pruning removes SNPs which are correlated 
# (measured by the squared correlation r2). From a group of correlated SNPs it retains the one with the 
# lowest $p$-value for the exposure. Use the function ieugwasr::ld_clump to prune the data. 
#The algorithm needs to know the rs identifier of the genetic variants (labeled as rsid) and 
#the $p$-value of the risk factor or exposure (labeled as pval).
IL33.merge_AD.clump = ieugwasr::ld_clump(IL33.merge_AD, clump_r2 = 0.01)
dim(IL33.merge_AD.clump)

ST2.merge_AD1.clump = ieugwasr::ld_clump(ST2.merge_AD1, clump_r2 = 0.01)
dim(ST2.merge_AD1.clump)

ST2.merge_AD2.clump = ieugwasr::ld_clump(ST2.merge_AD2, clump_r2 = 0.01)
dim(ST2.merge_AD2.clump)

# check why there is one SNP missing for ST2 with lower GWAS significance level
load("../GWAS_Data/iv_out.RData")
setdiff(ST2_clump_01$rsid, ST2.merge_AD1.clump$rsid)
'rs115540952' %in% AD.data$rsid # this SNP doesn't exist in AD GWAS so missed in ST2.merge_AD1.clump

#Calculate instrument strength
IL33.F=mean((IL33.merge_AD$beta_not_aligned/IL33.merge_AD$se)^2)
IL33.F
ST2.F1=mean((ST2.merge_AD1$beta_not_aligned/ST2.merge_AD1$se)^2)
ST2.F1
ST2.F2=mean((ST2.merge_AD2$beta_not_aligned/ST2.merge_AD2$se)^2)
ST2.F2

#################################################################################
# 5. Run the Mendelian Randomization 
#################################################################################
# IL33
IL33.rs = IL33.merge_AD.clump$rsid
IL33.beta = IL33.merge_AD.clump$beta_not_aligned
IL33.se = IL33.merge_AD.clump$se
IL33.AD_beta = IL33.merge_AD.clump$AD_beta
IL33.AD_se = IL33.merge_AD.clump$AD_se

IL33.mr.input = mr_input(bx = IL33.beta, bxse = IL33.se, by = IL33.AD_beta, byse = IL33.AD_se, 
                    exposure = "IL33", outcome = "Alzheimer's Disease", snps = IL33.rs)

mr_ivw(IL33.mr.input)

pdf(file = "../plots/MR_IVW_IL33.pdf")
mr_plot(IL33.mr.input, interactive=FALSE)
dev.off()

# the below line cannot run because: Method requires data on >2 variants
mr_allmethods(IL33.mr.input)
pdf(file = "../plots/MR_all_IL33.pdf")
mr_plot(mr_allmethods(IL33.mr.input, method = "main"))
dev.off()

# ST2 25IV
ST2.rs.1 = ST2.merge_AD1.clump$rsid
ST2.beta.1 = ST2.merge_AD1.clump$beta_not_aligned
ST2.se.1 = ST2.merge_AD1.clump$se
ST2.AD_beta.1 = ST2.merge_AD1.clump$AD_beta
ST2.AD_se.1 = ST2.merge_AD1.clump$AD_se

ST2.mr.input.1 = mr_input(bx = ST2.beta.1, bxse = ST2.se.1, by = ST2.AD_beta.1, byse = ST2.AD_se.1, 
                         exposure = "ST2-25IV", outcome = "Alzheimer's Disease", snps = ST2.rs.1)
mr_ivw(ST2.mr.input.1)

pdf(file = "../plots/MR_IVW_ST2_25IV.pdf")
mr_plot(ST2.mr.input.1, interactive=FALSE)
dev.off()

mr_allmethods(ST2.mr.input.1)
pdf(file = "../plots/MR_all_ST2_25IV.pdf")
mr_plot(mr_allmethods(ST2.mr.input.1, method = "main"))
dev.off()

# ST2 2IV
ST2.rs.2 = ST2.merge_AD2.clump$rsid
ST2.beta.2 = ST2.merge_AD2.clump$beta_not_aligned
ST2.se.2 = ST2.merge_AD2.clump$se
ST2.AD_beta.2 = ST2.merge_AD2.clump$AD_beta
ST2.AD_se.2 = ST2.merge_AD2.clump$AD_se

ST2.mr.input.2 = mr_input(bx = ST2.beta.2, bxse = ST2.se.2, by = ST2.AD_beta.2, byse = ST2.AD_se.2, 
                          exposure = "ST2-2IV", outcome = "Alzheimer's Disease", snps = ST2.rs.2)
mr_ivw(ST2.mr.input.2)

pdf(file = "../plots/MR_IVW_ST2_2IV.pdf")
mr_plot(ST2.mr.input.2, interactive=FALSE)
dev.off()

# the below line cannot run because: Method requires data on >2 variants
mr_allmethods(ST2.mr.input.2)
pdf(file = "../plots/MR_all_ST2_2IV.pdf")
mr_plot(mr_allmethods(ST2.mr.input.2, method = "main"))
dev.off()





