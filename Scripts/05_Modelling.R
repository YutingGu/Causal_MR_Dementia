# Load required library
library(dplyr)
library(parallel)
library(ggplot2)

#no_cores=detectCores()-1
#cl <- makeCluster(no_cores) 
#################################################################################
# Load Data
setwd("~/Desktop/MSc_HDAML_IC/course/SEM-2/TDS/Project/Rserver/Scripts")

data <- readRDS("../Data_Extraction/02_Cleaned_Data/Cleaned_ukb_with_PRS_new.rds")

Proteins_IL1RL1 <- readRDS("../Data_Extraction/02_Cleaned_Data/Proteins_IL1RL1_imputed.rds")
Proteins_IL33 <- readRDS("../Data_Extraction/02_Cleaned_Data/Proteins_IL33_imputed.rds")
Proteins_IL33_nd_filtered <- readRDS("../Data_Extraction/02_Cleaned_Data/Proteins_IL33_nd_filtered.rds")
Proteins_IL33_nd_filtered = Proteins_IL33_nd_filtered |> dplyr::rename(IL33_nd_filtered=IL33)
#################################################################################
# Combine all data together
data_full <- data
data_full <- left_join(data_full, Proteins_IL1RL1)
data_full <- left_join(data_full, Proteins_IL33)
data_full <- left_join(data_full, Proteins_IL33_nd_filtered)

#################################################################################
# Stratify dataset by gender
data_full_male <- subset(data_full, data_full$gender=='Male')
data_full_female <- subset(data_full, data_full$gender=='Female')

#################################################################################
# Define a function to extract parameters from a model - only for one exposure model
extract_parameters <- function(model){
  CI <- exp(confint(model))
  model_summary = summary(model)
  
  model_data <- as.character(model_summary[["call"]][["data"]])
  outcome <- as.character(model_summary[["call"]][["formula"]][[2]])
  exposure <- as.character(model_summary[["call"]][["formula"]][[3]])
  OR <- unname(exp(model$coefficients)[2])
  lower <- CI[2,1]
  upper <- CI[2,2]
  pvalue <- model_summary[["coefficients"]][2,4]
  
  c(model_data, outcome, exposure, OR, lower, upper, pvalue)
}

extract_parameters2 <- function(model){
  CI <- exp(confint(model))
  model_summary = summary(model)
  
  model_data <- as.character(model_summary[["call"]][["data"]])
  outcome <- as.character(model_summary[["call"]][["formula"]][[2]])
  exposure1 <- as.character(model_summary[["call"]][["formula"]][[3]][[2]])
  exposure2 <- as.character(model_summary[["call"]][["formula"]][[3]][[3]])
  exposures = paste(exposure1,'+',exposure2)
  
  OR1 <- unname(exp(model$coefficients)[2])
  OR2 <- unname(exp(model$coefficients)[3])
  lower1 <- CI[2,1]
  upper1 <- CI[2,2]
  lower2 <- CI[3,1]
  upper2 <- CI[3,2]
  pvalue1 <- model_summary[["coefficients"]][2,4]
  pvalue2 <- model_summary[["coefficients"]][3,4]
  
  
  variable1 <- c(model_data, outcome, exposures, exposure1, OR1, lower1, upper1, pvalue1)
  variable2 <- c(model_data, outcome, exposures, exposure2, OR2, lower2, upper2, pvalue2)
  rbind(variable1,variable2)
}


#################################################################################
# Some notes of the models:
## 1. glm will automatic ignore rows with outcome being NA, so it doesn't matter if outcome contains NA
#################################################################################
# Modelling - Observational study: Dementia ~ Protein

column_names <- c("model_data", "outcome", "exposure", "OR", "lower", "upper", "pvalue")
ST2_obs_results <- data.frame(matrix(ncol = length(column_names), nrow = 0))

## ST2
Obs.model_all.ST2 <- glm(case_all ~ IL1RL1, data=data_full, family = binomial(link = 'logit'))
ST2_obs_results <- rbind(ST2_obs_results, extract_parameters(Obs.model_all.ST2))
Obs.model_AD.ST2 <- glm(case_AD ~ IL1RL1, data=data_full, family = binomial(link = 'logit'))
ST2_obs_results <- rbind(ST2_obs_results, extract_parameters(Obs.model_AD.ST2))
Obs.model_VD.ST2 <- glm(case_VD ~ IL1RL1, data=data_full, family = binomial(link = 'logit'))
ST2_obs_results <- rbind(ST2_obs_results, extract_parameters(Obs.model_VD.ST2))
Obs.model_OD.ST2 <- glm(case_OD ~ IL1RL1, data=data_full, family = binomial(link = 'logit'))
ST2_obs_results <- rbind(ST2_obs_results, extract_parameters(Obs.model_OD.ST2))

colnames(ST2_obs_results) <- column_names
write.csv(ST2_obs_results, "../results/ST2_obs_results.csv", row.names = FALSE)

plot(data_full$case_all, data_full$IL1RL1)
boxplot(data_full[data_full$case_all==0,]$IL1RL1, data_full[data_full$case_all==1,]$IL1RL1)
plot(data_full$case_AD, data_full$IL1RL1)
boxplot(data_full[data_full$case_AD==0,]$IL1RL1, data_full[data_full$case_AD==1,]$IL1RL1)
plot(data_full$case_VD, data_full$IL1RL1)
boxplot(data_full[data_full$case_VD==0,]$IL1RL1, data_full[data_full$case_VD==1,]$IL1RL1)
plot(data_full$case_OD, data_full$IL1RL1)
boxplot(data_full[data_full$case_OD==0,]$IL1RL1, data_full[data_full$case_OD==1,]$IL1RL1)


## IL33
column_names <- c("model_data", "outcome", "exposure", "OR", "lower", "upper", "pvalue")
IL33_obs_results <- data.frame(matrix(ncol = length(column_names), nrow = 0))

Obs.model_all.IL33 <- glm(case_all ~ IL33, data=data_full, family = binomial(link = 'logit'))
IL33_obs_results <- rbind(IL33_obs_results, extract_parameters(Obs.model_all.IL33))
Obs.model_AD.IL33 <- glm(case_AD ~ IL33, data=data_full, family = binomial(link = 'logit'))
IL33_obs_results <- rbind(IL33_obs_results, extract_parameters(Obs.model_AD.IL33))
Obs.model_VD.IL33 <- glm(case_VD ~ IL33, data=data_full, family = binomial(link = 'logit'))
IL33_obs_results <- rbind(IL33_obs_results, extract_parameters(Obs.model_VD.IL33))
Obs.model_OD.IL33 <- glm(case_OD ~ IL33, data=data_full, family = binomial(link = 'logit'))
IL33_obs_results <- rbind(IL33_obs_results, extract_parameters(Obs.model_OD.IL33))

colnames(IL33_obs_results) <- column_names
write.csv(IL33_obs_results, "../results/IL33_obs_results.csv", row.names = FALSE)

plot(data_full$case_all, data_full$IL33)
boxplot(data_full[data_full$case_all==0,]$IL33, data_full[data_full$case_all==1,]$IL33)
plot(data_full$case_AD, data_full$IL33)
boxplot(data_full[data_full$case_AD==0,]$IL33, data_full[data_full$case_AD==1,]$IL33)
plot(data_full$case_VD, data_full$IL33)
boxplot(data_full[data_full$case_VD==0,]$IL33, data_full[data_full$case_VD==1,]$IL33)
plot(data_full$case_OD, data_full$IL33)
boxplot(data_full[data_full$case_OD==0,]$IL33, data_full[data_full$case_OD==1,]$IL33)

#################################################################################
# Modelling - PRS checking: Risk factor ~ PRS
PRS_RF.model.PRS_ST2_26IV <- lm(data_full$IL1RL1 ~ data_full$PRS_ST2_26IV)
summary(PRS_RF.model.PRS_ST2_26IV)
PRS_RF.model.PRS_ST2_2IV <- lm(data_full$IL1RL1 ~ data_full$PRS_ST2_2IV)
summary(PRS_RF.model.PRS_ST2_2IV)
PRS_RF.model.PRS_IL33_2IV <- lm(data_full$IL33 ~ data_full$PRS_IL33_2IV)
summary(PRS_RF.model.PRS_IL33_2IV)
PRS_RF.model.PRS_IL33_2IV_nd <- lm(data_full$IL33_nd_filtered ~ data_full$PRS_IL33_2IV)
summary(PRS_RF.model.PRS_IL33_2IV_nd)

plot(data_full$PRS_ST2_26IV, data_full$IL1RL1)
plot(data_full$PRS_ST2_2IV, data_full$IL1RL1)
plot(data_full$PRS_IL33_2IV, data_full$IL33)
plot(data_full$PRS_IL33_2IV, data_full$IL33_nd_filtered)
#################################################################################
# Modelling - Positive Control: Asthma ~ PRS
# number of cases
table(data_full$case_asthma)

table(data_full$case_asthma, useNA = 'ifany')
PC.model.PRS_ST2_26IV <- glm(case_asthma ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
summary(PC.model.PRS_ST2_26IV)
exp(PC.model.PRS_ST2_26IV$coefficients)
exp(confint(PC.model.PRS_ST2_26IV))

PC.model.PRS_ST2_2IV <- glm(case_asthma ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
summary(PC.model.PRS_ST2_2IV)
exp(PC.model.PRS_ST2_2IV$coefficients)
exp(confint(PC.model.PRS_ST2_2IV))

PC.model.PRS_IL33_2IV <- glm(case_asthma ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
summary(PC.model.PRS_IL33_2IV)
exp(PC.model.PRS_IL33_2IV$coefficients)
exp(confint(PC.model.PRS_IL33_2IV))

#################################################################################
# Modelling - Negative Control: HotDrinkTemp_Binary ~ PRS
data_full[,c('PRS_ST2_26IV','HotDrinkTemp_Binary')] |> 
  ggplot(aes(PRS_ST2_26IV, color = HotDrinkTemp_Binary)) + 
  geom_density() + 
  labs(title = "Distribution of PRS", 
       subtitle = "Outcome = HotDrinkTemp_Binary\tPRS = ST2 - 26 SNPs", 
       x="PRS", color = "") 

table(data_full$HotDrinkTemp_Binary, useNA = 'ifany')

NC.model.PRS_ST2_26IV <- glm(HotDrinkTemp_Binary ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
summary(NC.model.PRS_ST2_26IV)
exp(NC.model.PRS_ST2_26IV$coefficients)
exp(confint(NC.model.PRS_ST2_26IV))

NC.model.PRS_ST2_2IV <- glm(HotDrinkTemp_Binary ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
summary(NC.model.PRS_ST2_2IV)
exp(NC.model.PRS_ST2_2IV$coefficients)
exp(confint(NC.model.PRS_ST2_2IV))

NC.model.PRS_IL33_2IV <- glm(HotDrinkTemp_Binary ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
summary(NC.model.PRS_IL33_2IV)
exp(NC.model.PRS_IL33_2IV$coefficients)
exp(confint(NC.model.PRS_IL33_2IV))

#################################################################################
# Modelling - MR: Dementia ~ PRS - No stratification
# number of cases
table(data_full$case_all)
table(data_full$case_AD)
table(data_full$case_VD)
table(data_full$case_OD)

# Define a empty dataframe to store results
column_names <- c("model_data", "outcome", "exposure", "OR", "lower", "upper", "pvalue")
MR_results_nostrat <- data.frame(matrix(ncol = length(column_names), nrow = 0))

# All
## Logistic: case_all ~ PRS_ST2_26IV
MR.model_all.PRS_ST2_26IV <- glm(case_all ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_all.PRS_ST2_26IV))
## Logistic: case_all ~ PRS_ST2_2IV
MR.model_all.PRS_ST2_2IV <- glm(case_all ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_all.PRS_ST2_2IV))
## Logistic: case_all ~ PRS_IL33_2IV
MR.model_all.PRS_IL33_2IV <- glm(case_all ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_all.PRS_IL33_2IV))

# AD
## Logistic: case_AD ~ PRS_ST2_26IV
MR.model_AD.PRS_ST2_26IV <- glm(case_AD ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_AD.PRS_ST2_26IV))
## Logistic: case_AD ~ PRS_ST2_2IV
MR.model_AD.PRS_ST2_2IV <- glm(case_AD ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_AD.PRS_ST2_2IV))
## Logistic: case_AD ~ PRS_IL33_2IV
MR.model_AD.PRS_IL33_2IV <- glm(case_AD ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_AD.PRS_IL33_2IV))

# VD
## Logistic: case_VD ~ PRS_ST2_26IV
MR.model_VD.PRS_ST2_26IV <- glm(case_VD ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_VD.PRS_ST2_26IV))
## Logistic: case_VD ~ PRS_ST2_2IV
MR.model_VD.PRS_ST2_2IV <- glm(case_VD ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_VD.PRS_ST2_2IV))
## Logistic: case_VD ~ PRS_IL33_2IV
MR.model_VD.PRS_IL33_2IV <- glm(case_VD ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_VD.PRS_IL33_2IV))

# OD
## Logistic: case_OD ~ PRS_ST2_26IV
MR.model_OD.PRS_ST2_26IV <- glm(case_OD ~ PRS_ST2_26IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_OD.PRS_ST2_26IV))
## Logistic: case_OD ~ PRS_ST2_2IV
MR.model_OD.PRS_ST2_2IV <- glm(case_OD ~ PRS_ST2_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_OD.PRS_ST2_2IV))
## Logistic: case_OD ~ PRS_IL33_2IV
MR.model_OD.PRS_IL33_2IV <- glm(case_OD ~ PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_nostrat <- rbind(MR_results_nostrat, extract_parameters(MR.model_OD.PRS_IL33_2IV))

colnames(MR_results_nostrat) <- column_names
write.csv(MR_results_nostrat, "../results/MR_results_nostrat", row.names = FALSE)
#################################################################################
# Modelling - MR: Dementia ~ PRS_ST2 + PRS_IL33 - No stratification

# Define a empty dataframe to store results
column_names <- c("model_data", "outcome", "exposures", "exposure", "OR", "lower", "upper", "pvalue")
MR_results_conditioned <- data.frame(matrix(ncol = length(column_names), nrow = 0))

# All
## Logistic: case_all ~ PRS_ST2_26IV + PRS_IL33
MR.model_all.PRS_ST2_26IV_PRS_IL33_2IV <- glm(case_all ~ PRS_ST2_26IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_all.PRS_ST2_26IV_PRS_IL33_2IV))
## Logistic: case_all ~ PRS_ST2_2IV + PRS_IL33
MR.model_all.PRS_ST2_2IV_PRS_IL33_2IV <- glm(case_all ~ PRS_ST2_2IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_all.PRS_ST2_2IV_PRS_IL33_2IV))

# AD
## Logistic: case_AD ~ PRS_ST2_26IV + PRS_IL33
MR.model_AD.PRS_ST2_26IV_PRS_IL33_2IV <- glm(case_AD ~ PRS_ST2_26IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_AD.PRS_ST2_26IV_PRS_IL33_2IV))
## Logistic: case_AD ~ PRS_ST2_2IV + PRS_IL33
MR.model_AD.PRS_ST2_2IV_PRS_IL33_2IV <- glm(case_AD ~ PRS_ST2_2IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_AD.PRS_ST2_2IV_PRS_IL33_2IV))

# AD
## Logistic: case_VD ~ PRS_ST2_26IV + PRS_IL33
MR.model_VD.PRS_ST2_26IV_PRS_IL33_2IV <- glm(case_VD ~ PRS_ST2_26IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_VD.PRS_ST2_26IV_PRS_IL33_2IV))
## Logistic: case_VD ~ PRS_ST2_2IV + PRS_IL33
MR.model_VD.PRS_ST2_2IV_PRS_IL33_2IV <- glm(case_VD ~ PRS_ST2_2IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_VD.PRS_ST2_2IV_PRS_IL33_2IV))

# OD
## Logistic: case_OD ~ PRS_ST2_26IV + PRS_IL33
MR.model_OD.PRS_ST2_26IV_PRS_IL33_2IV <- glm(case_OD ~ PRS_ST2_26IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_OD.PRS_ST2_26IV_PRS_IL33_2IV))
## Logistic: case_OD ~ PRS_ST2_2IV + PRS_IL33
MR.model_OD.PRS_ST2_2IV_PRS_IL33_2IV <- glm(case_OD ~ PRS_ST2_2IV + PRS_IL33_2IV, data=data_full, family = binomial(link = 'logit'))
MR_results_conditioned <- rbind(MR_results_conditioned, extract_parameters2(MR.model_OD.PRS_ST2_2IV_PRS_IL33_2IV))

colnames(MR_results_conditioned) <- column_names
write.csv(MR_results_conditioned, "../results/MR_results_conditioned.csv", row.names = FALSE)

#################################################################################
# Modelling - MR: Dementia ~ PRS - Sex stratification
# number of cases
table(data_full_male$case_all)
table(data_full_female$case_all)
table(data_full_male$case_AD)
table(data_full_female$case_AD)
table(data_full_male$case_VD)
table(data_full_female$case_VD)
table(data_full_male$case_OD)
table(data_full_female$case_OD)

# Define a empty dataframe to store results
column_names <- c("model_data", "outcome", "exposure", "OR", "lower", "upper", "pvalue")
MR_results_gender <- data.frame(matrix(ncol = length(column_names), nrow = 0))

# All
## Logistic: case_all ~ PRS_ST2_26IV
MR_Male.model_all.PRS_ST2_26IV <- glm(case_all ~ PRS_ST2_26IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_all.PRS_ST2_26IV))

MR_Female.model_all.PRS_ST2_26IV <- glm(case_all ~ PRS_ST2_26IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_all.PRS_ST2_26IV))

## Logistic: case_all ~ PRS_ST2_2IV
MR_Male.model_all.PRS_ST2_2IV <- glm(case_all ~ PRS_ST2_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_all.PRS_ST2_2IV))

MR_Female.model_all.PRS_ST2_2IV <- glm(case_all ~ PRS_ST2_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_all.PRS_ST2_2IV))

## Logistic: case_all ~ PRS_IL33_2IV
MR_Male.model_all.PRS_IL33_2IV <- glm(case_all ~ PRS_IL33_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_all.PRS_IL33_2IV))

MR_Female.model_all.PRS_IL33_2IV <- glm(case_all ~ PRS_IL33_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_all.PRS_IL33_2IV))


# AD
## Logistic: case_AD ~ PRS_ST2_26IV
MR_Male.model_AD.PRS_ST2_26IV <- glm(case_AD ~ PRS_ST2_26IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_AD.PRS_ST2_26IV))

MR_Female.model_AD.PRS_ST2_26IV <- glm(case_AD ~ PRS_ST2_26IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_AD.PRS_ST2_26IV))

## Logistic: case_AD ~ PRS_ST2_2IV
MR_Male.model_AD.PRS_ST2_2IV <- glm(case_AD ~ PRS_ST2_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_AD.PRS_ST2_2IV))

MR_Female.model_AD.PRS_ST2_2IV <- glm(case_AD ~ PRS_ST2_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_AD.PRS_ST2_2IV))

## Logistic: case_AD ~ PRS_IL33_2IV
MR_Male.model_AD.PRS_IL33_2IV <- glm(case_AD ~ PRS_IL33_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_AD.PRS_IL33_2IV))

MR_Female.model_AD.PRS_IL33_2IV <- glm(case_AD ~ PRS_IL33_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_AD.PRS_IL33_2IV))


# VD
## Logistic: case_VD ~ PRS_ST2_26IV
MR_Male.model_VD.PRS_ST2_26IV <- glm(case_VD ~ PRS_ST2_26IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_VD.PRS_ST2_26IV))

MR_Female.model_VD.PRS_ST2_26IV <- glm(case_VD ~ PRS_ST2_26IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_VD.PRS_ST2_26IV))

## Logistic: case_VD ~ PRS_ST2_2IV
MR_Male.model_VD.PRS_ST2_2IV <- glm(case_VD ~ PRS_ST2_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_VD.PRS_ST2_2IV))

MR_Female.model_VD.PRS_ST2_2IV <- glm(case_VD ~ PRS_ST2_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_VD.PRS_ST2_2IV))

## Logistic: case_VD ~ PRS_IL33_2IV
MR_Male.model_VD.PRS_IL33_2IV <- glm(case_VD ~ PRS_IL33_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_VD.PRS_IL33_2IV))

MR_Female.model_VD.PRS_IL33_2IV <- glm(case_VD ~ PRS_IL33_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_VD.PRS_IL33_2IV))


# OD
## Logistic: case_OD ~ PRS_ST2_26IV
MR_Male.model_OD.PRS_ST2_26IV <- glm(case_OD ~ PRS_ST2_26IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_OD.PRS_ST2_26IV))

MR_Female.model_OD.PRS_ST2_26IV <- glm(case_OD ~ PRS_ST2_26IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_OD.PRS_ST2_26IV))

## Logistic: case_OD ~ PRS_ST2_2IV
MR_Male.model_OD.PRS_ST2_2IV <- glm(case_OD ~ PRS_ST2_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_OD.PRS_ST2_2IV))

MR_Female.model_OD.PRS_ST2_2IV <- glm(case_OD ~ PRS_ST2_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_OD.PRS_ST2_2IV))

## Logistic: case_OD ~ PRS_IL33_2IV
MR_Male.model_OD.PRS_IL33_2IV <- glm(case_OD ~ PRS_IL33_2IV, data=data_full_male, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Male.model_OD.PRS_IL33_2IV))

MR_Female.model_OD.PRS_IL33_2IV <- glm(case_OD ~ PRS_IL33_2IV, data=data_full_female, family = binomial(link = 'logit'))
MR_results_gender <- rbind(MR_results_gender, extract_parameters(MR_Female.model_OD.PRS_IL33_2IV))

colnames(MR_results_gender) <- column_names
write.csv(MR_results_gender, "../results/MR_results_gender.csv", row.names = FALSE)
