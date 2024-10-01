# Load required library
library(dplyr)
library(parallel)
library(table1)
library(ggplot2)

#no_cores=detectCores()-1
#cl <- makeCluster(no_cores) 
#################################################################################
# Load Data
setwd("/rds/general/project/hda_23-24/live/TDS/group1/Scripts")

## load processed ukb dataset with lables
data_full <- readRDS("../Data_Extraction/02_Cleaned_Data/Cleaned_ukb_with_outcome.rds")

## genetic data ST2
genetic_data_ST2 <- readRDS("../Data_Extraction/02_Cleaned_Data/genetic_data_numeric_ST2.rds")
genetic_data_ST2 <- as.data.frame(apply(genetic_data_ST2, c(1,2), as.integer))
genetic_data_ST2$eid_ST2 <- rownames(genetic_data_ST2)

## genetic data IL33
genetic_data_IL33 <- readRDS("../Data_Extraction/02_Cleaned_Data/genetic_data_numeric_IL33.rds")
genetic_data_IL33 <- as.data.frame(apply(genetic_data_IL33, c(1,2), as.integer))
genetic_data_IL33$eid_IL33 <- rownames(genetic_data_IL33)

## load IV GWAS
IL33_GWAS <- readRDS("../GWAS_Data/snp_table_group1_IL33.rds")
ST2_GWAS <- readRDS("../GWAS_Data/snp_table_group1_ST2.rds")
#################################################################################
#genetic_data_ST2 = genetic_data_ST2 |> dplyr::rename(eid_ST2=eid)
#genetic_data_IL33 = genetic_data_IL33 |> dplyr::rename(eid_IL33=eid)

# join genotype data
data_full <- left_join(data_full, genetic_data_ST2, join_by(eid==eid_ST2), keep=TRUE)
data_full <- left_join(data_full, genetic_data_IL33, join_by(eid==eid_IL33), keep=TRUE)

# check how many data are not matched
data_full$genotype_NA <- ifelse(rowSums(is.na(data_full[ST2_SNP_list])) == 26 | rowSums(is.na(data_full[IL33_SNP_list])) == 2, TRUE, FALSE)
table(data_full$genotype_NA)


#################################################################################
# IL33_beta <- IL33_clump_01 |> select(rsid,beta)
# ST2_beta <- ST2_clump_01 |> select(rsid,beta)

IL33_beta <- IL33_GWAS |> select(rsid,beta_not_aligned)
ST2_beta <- ST2_GWAS |> select(rsid,beta_not_aligned)

# create SNP-id list for each protein
ST2_SNP_list <- names(genetic_data_ST2)[1:length(names(genetic_data_ST2))-1]
ST2_SNP_list_short <- c('rs13029918','rs1468789')
IL33_SNP_list <- names(genetic_data_IL33)[1:length(names(genetic_data_IL33))-1]

# check how many participants having missing genotype data for each SNP 
colSums(data_full[ST2_SNP_list]==0 | is.na(data_full[ST2_SNP_list]))
colSums(data_full[ST2_SNP_list_short]==0 | is.na(data_full[ST2_SNP_list_short]))
colSums(data_full[IL33_SNP_list]==0 | is.na(data_full[IL33_SNP_list]))

# compute the number of individuals of each number of missing variants
table(rowSums(is.na(data_full[ST2_SNP_list])))
table(rowSums(is.na(data_full[ST2_SNP_list_short])))
table(rowSums(is.na(data_full[IL33_SNP_list])))

# exclude all NA rows
data_genotype_NA <- data_full[data_full$genotype_NA,]
data_genotype_NA <- subset(data_genotype_NA, select = -c(eid_ST2,eid_IL33,genotype_NA))
data_full <- data_full[!data_full$genotype_NA,]
data_full <- subset(data_full, select = -c(eid_ST2,eid_IL33,genotype_NA))

# make ST2_beta dataset in a format which is easy to compute
rownames(ST2_beta) <- ST2_beta[, 1]
ST2_beta <- ST2_beta[ST2_SNP_list,]

# make IL33_beta dataset in a format which is easy to compute
IL33_beta <- as.data.frame(IL33_beta)
rownames(IL33_beta) <- IL33_beta[, 1]
IL33_beta <- IL33_beta[IL33_SNP_list,]

# compute number of variants to divide when computing PRS
number_variants_ST2_26IV <- 26 - rowSums(is.na(data_full[ST2_SNP_list]))
number_variants_ST2_2IV <- 2 - rowSums(is.na(data_full[ST2_SNP_list_short]))
number_variants_IL33_2IV <- 2 - rowSums(is.na(data_full[IL33_SNP_list]))
## check data
table(number_variants_ST2_26IV, useNA = 'ifany')
table(number_variants_ST2_2IV, useNA = 'ifany')
table(number_variants_IL33_2IV, useNA = 'ifany')

# convert NA genotype to 0
apply(data_full[,c(ST2_SNP_list,IL33_SNP_list)], 2, FUN = function(x) {table(x, useNA = 'ifany')})
data_full[,c(ST2_SNP_list)] <- replace(data_full[,c(ST2_SNP_list)],is.na(data_full[,c(ST2_SNP_list)]),0)
data_full[,c(IL33_SNP_list)] <- replace(data_full[,c(IL33_SNP_list)],is.na(data_full[,c(IL33_SNP_list)]),0)

## check data
apply(data_full[,c(ST2_SNP_list,IL33_SNP_list)], 2, FUN = function(x) {table(x, useNA = 'ifany')})

# Compute PRS
## ST2 26 variants
data_full$PRS_ST2_26IV <- (as.matrix(data_full[,ST2_SNP_list]) %*% ST2_beta[,2]) / number_variants_ST2_26IV
## ST2 2 variants
data_full$PRS_ST2_2IV <- (as.matrix(data_full[,ST2_SNP_list_short]) %*% ST2_beta[ST2_SNP_list_short,2]) / number_variants_ST2_2IV
## IL33 2 variants
data_full$PRS_IL33_2IV <- (as.matrix(data_full[,IL33_SNP_list]) %*% IL33_beta[,2]) / number_variants_IL33_2IV
data_full$PRS_IL33_2IV <- replace(data_full$PRS_IL33_2IV,is.na(data_full$PRS_IL33_2IV),0) # due to zero division, need to correct computing error
table(data_full$PRS_IL33_2IV, useNA = 'ifany')

#################################################################################
# Cross check with Lucie's result
#Lucie_data_ST2_26_outcomes_asthma_and_PRS <- readRDS("/rds/general/project/hda_23-24/live/TDS/group1/Data_Extraction/02_Cleaned_Data/Lucie_data_ST2_26_outcomes_asthma_and_PRS.rds")

#PRS_check <- left_join(data_full[c('eid','PRS_ST2_26IV')], Lucie_data_ST2_26_outcomes_asthma_and_PRS[c('eid','PRS_ST2_26')])
#PRS_check$diff = PRS_check$PRS_ST2_26IV - PRS_check$PRS_ST2_26
#table(PRS_check$diff, useNA = 'ifany')

#PRS_check$diff_binary <- ifelse(abs(PRS_check$diff) <= 0.00001,TRUE,FALSE)
#table(PRS_check$diff_binary, useNA = 'ifany')

#################################################################################
# Save data
saveRDS(data_full, "../Data_Extraction/02_Cleaned_Data/Cleaned_ukb_with_PRS_new.rds")

#################################################################################
# Missingness plots
# compute the number of individuals of each number of missing variants
PRS_missing_distribution.plot <- function (tab,title) {
  print("raw numbers:")
  print(cbind(tab,rowSums(tab)))
  
  table_p = round(tab/rowSums(tab) * 100,2)
  print("Percentage:")
  print(table_p)
  
  table_p = as.data.frame(table_p)
  colnames(table_p) = c("case","number","percent")
  table_p |> 
    mutate(case=ifelse(case==TRUE,"Case","Control")) |>
    ggplot(aes(x = number, y=percent, fill = case)) + 
    geom_col(position = "dodge") + 
    labs(title = "Genotype Missingness Distribution", 
         subtitle = title,
         x="Number of missing variants", 
         y="Population Percentage", 
         fill = "") 
}



# all types
tab <- table(as.logical(data_full$case_all),26 - number_variants_ST2_26IV)
PRS_missing_distribution.plot(tab,"Outcome = All types\tPRS = ST2 - 26 SNPs")
tab <- table(as.logical(data_full$case_all),2 - number_variants_ST2_2IV)
PRS_missing_distribution.plot(tab,"Outcome = All types\tPRS = ST2 - 2 SNPs")
tab <- table(as.logical(data_full$case_all),2 - number_variants_IL33_2IV)
PRS_missing_distribution.plot(tab,"Outcome = All types\tPRS = IL33 - 2 SNPs")

# AD
tab <- table(as.logical(data_full$case_AD),26 - number_variants_ST2_26IV)
PRS_missing_distribution.plot(tab,"Outcome = Alzheimer's Disease\tPRS = ST2 - 26 SNPs")
tab <- table(as.logical(data_full$case_AD),2 - number_variants_ST2_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Alzheimer's Disease\tPRS = ST2 - 2 SNPs")
tab <- table(as.logical(data_full$case_AD),2 - number_variants_IL33_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Alzheimer's Disease\tPRS = IL33 - 2 SNPs")

# VD
tab <- table(as.logical(data_full$case_VD),26 - number_variants_ST2_26IV)
PRS_missing_distribution.plot(tab,"Outcome = Vascular Dementia\tPRS = ST2 - 26 SNPs")
tab <- table(as.logical(data_full$case_VD),2 - number_variants_ST2_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Vascular Dementia\tPRS = ST2 - 2 SNPs")
tab <- table(as.logical(data_full$case_VD),2 - number_variants_IL33_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Vascular Dementia\tPRS = IL33 - 2 SNPs")

# OD
tab <- table(as.logical(data_full$case_OD),26 - number_variants_ST2_26IV)
PRS_missing_distribution.plot(tab,"Outcome = Other Dementia\tPRS = ST2 - 26 SNPs")
tab <- table(as.logical(data_full$case_OD),2 - number_variants_ST2_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Other Dementia\tPRS = ST2 - 2 SNPs")
tab <- table(as.logical(data_full$case_OD),2 - number_variants_IL33_2IV)
PRS_missing_distribution.plot(tab,"Outcome = Other Dementia\tPRS = IL33 - 2 SNPs")


#################################################################################
# Visualisation

## plot PRS
plot_data <- data_full[c('PRS_ST2_26IV','PRS_ST2_2IV','PRS_IL33_2IV','case_all', 'case_AD','case_VD','case_OD')]
plot_data <- plot_data[complete.cases(plot_data),]

PRS_density.plot <- function(data,title) {
  colnames(data) = c("PRS","case")
  data |> 
    mutate(case=ifelse(case==TRUE|case==1,"Case","Control")) |>
    ggplot(aes(PRS, color = case)) + 
    geom_density() + 
    labs(title = "Distribution of PRS", subtitle = title, x="PRS", color = "") 
}


# all types
PRS_density.plot(plot_data[,c('PRS_ST2_26IV','case_all')], "Outcome = All types\tPRS = ST2 - 26 SNPs")
PRS_density.plot(plot_data[,c('PRS_ST2_2IV','case_all')], "Outcome = All types\tPRS = ST2 - 2 SNPs")
PRS_density.plot(plot_data[,c('PRS_IL33_2IV','case_all')], "Outcome = All types\tPRS = IL33 - 2 SNPs")

# AD
PRS_density.plot(plot_data[,c('PRS_ST2_26IV','case_AD')], "Outcome = Alzheimer's Disease\tPRS = ST2 - 26 SNPs")
PRS_density.plot(plot_data[,c('PRS_ST2_2IV','case_AD')], "Outcome = Alzheimer's Disease\tPRS = ST2 - 2 SNPs")
PRS_density.plot(plot_data[,c('PRS_IL33_2IV','case_AD')], "Outcome = Alzheimer's Disease\tPRS = IL33 - 2 SNPs")

# VD
PRS_density.plot(plot_data[,c('PRS_ST2_26IV','case_VD')], "Outcome = Vascular Dementia\tPRS = ST2 - 26 SNPs")
PRS_density.plot(plot_data[,c('PRS_ST2_2IV','case_VD')], "Outcome = Vascular Dementia\tPRS = ST2 - 2 SNPs")
PRS_density.plot(plot_data[,c('PRS_IL33_2IV','case_VD')], "Outcome = Vascular Dementia\tPRS = IL33 - 2 SNPs")

# Other
PRS_density.plot(plot_data[,c('PRS_ST2_26IV','case_OD')], "Outcome = Other Dementia\tPRS = ST2 - 26 SNPs")
PRS_density.plot(plot_data[,c('PRS_ST2_2IV','case_OD')], "Outcome = Other Dementia\tPRS = ST2 - 2 SNPs")
PRS_density.plot(plot_data[,c('PRS_IL33_2IV','case_OD')], "Outcome = Other Dementia\tPRS = IL33 - 2 SNPs")



#plot(density(x[x$case_all==TRUE,PRS], na.rm = TRUE), lwd = 2, col = "red", main = "", las = 1, xlab = "", ylab = "")
#lines(density(x[x$case_all==FALSE,PRS], na.rm = TRUE), lwd = 2, col = "blue")