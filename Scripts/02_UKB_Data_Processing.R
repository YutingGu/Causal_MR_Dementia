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
## Outcome
output_final_AD <- readRDS("../Data_Extraction/01_Extracted_Data/output_final_AD.rds")
output_final_VD <- readRDS("../Data_Extraction/01_Extracted_Data/output_final_VD.rds")
output_final_OD <- readRDS("../Data_Extraction/01_Extracted_Data/output_final_OD.rds")
output_final_self_report_dementia <- readRDS("../Data_Extraction/01_Extracted_Data/output_final_self_report_dementia.rds")

## Features
ukb_recoded <- readRDS("../Data_Extraction/01_Extracted_Data/ukb_recoded.rds")
ukb_recoded$eid <- rownames(ukb_recoded)
annot <- readRDS("../Data_Extraction/extraction_and_recoding/outputs/annot.rds")

# Positive and Negative Control
NegativeControl <- readRDS("../Data_Extraction/01_Extracted_Data/NegativeControl_HotDrink.rds")
output_final_asthma <- readRDS("../Data_Extraction/01_Extracted_Data/output_final_asthma.rds")

# try to compute in parallel
# genetic_data_ST2_test <- as.data.frame(parApply(cl = cl, X = genetic_data_ST2,MARGIN=c(1,2), FUN = as.integer))

#################################################################################
# Process columns
## 
ukb_processed <- ukb_recoded |> 
  mutate(
    LFU_flag = ifelse(is.na(`00_LostFollowUp_date.0.0`),0,1), # whether lost follow-up: 1 - lost, 0 - not lost,
    birth_year = `01_Birth_yr.0.0`,
    age_now = 2024 - `01_Birth_yr.0.0`,
    gender = `01_Sex.0.0`
    #age_diag = as.numeric(substr(date_diagnosis,1,4)) - `01_Birth_yr.0.0`,
    #age_recr = as.numeric(substr(date_recr,1,4)) - `01_Birth_yr.0.0`,
  ) |> 
  select(eid, LFU_flag,birth_year, age_now, gender)

## Process education level seperately
ukb_processed$edu_level <- ukb_recoded$`01_Edu_level.0.0` == 'College or University degree'
for (i in 1:length(ukb_processed)[1]) {
  edu_1 = ukb_recoded$`01_Edu_level.1.0`[i] == 'College or University degree'
  edu_2 = ukb_recoded$`01_Edu_level.2.0`[i] == 'College or University degree'
  edu_3 = ukb_recoded$`01_Edu_level.3.0`[i] == 'College or University degree'
  if (isTRUE(edu_1)) {
    ukb_processed$edu_level[i] = TRUE
  }else if (isTRUE(edu_2)){
    ukb_processed$edu_level[i] = TRUE
  }else if (isTRUE(edu_3)){
    ukb_processed$edu_level[i] = TRUE
  }
}
table(ukb_processed$edu_level, useNA = 'ifany')

## ethinicity group
ethnicity_list = c('White', 'British', 'Irish', 'Any other white background')
ukb_processed$ethnicity_white = ukb_recoded$`01_Ethnicity.0.0` %in% ethnicity_list 
for (i in 1:length(ukb_processed)[1]) {
  eth_1 = ukb_recoded$`01_Ethnicity.1.0`[i] %in% ethnicity_list
  eth_2 = ukb_recoded$`01_Ethnicity.2.0`[i] %in% ethnicity_list
  eth_3 = ukb_recoded$`01_Ethnicity.3.0`[i] %in% ethnicity_list
  if (isTRUE(eth_1)) {
    ukb_processed$ethnicity_white[i] = TRUE
  }else if (isTRUE(eth_2)){
    ukb_processed$ethnicity_white[i] = TRUE
  }else if (isTRUE(eth_3)){
    ukb_processed$ethnicity_white[i] = TRUE
  }
}
table(ukb_processed$ethnicity_white, useNA = 'ifany')

## SBP only use first instance, average two readings
ukb_processed$SBP = apply(ukb_recoded[,c('03_SBP.0.0','03_SBP.0.1')], 1, FUN = function(x) mean(x,na.rm = TRUE))
table(is.na(ukb_processed$`SBP` ), useNA = 'ifany')

## diabetes
# table(ukb_recoded$`04_Diabetes_doc.0.0`, useNA = 'ifany')
# ukb_processed$diabetes = ukb_recoded$`04_Diabetes_doc.0.0` == 'Yes'
  
## Total Cholesterol
table(is.na(ukb_recoded$`03_TotalChol.0.0`), useNA = 'ifany')
table(is.na(ukb_recoded$`03_TotalChol.1.0`), useNA = 'ifany')
ukb_processed$TotalChol = apply(ukb_recoded[,c('03_TotalChol.0.0','03_TotalChol.1.0')], 1, FUN = function(x) mean(x,na.rm = TRUE))
table(is.na(ukb_processed$`TotalChol` ), useNA = 'ifany')

## BMI
table(is.na(ukb_recoded$`03_BMI.0.0` ), useNA = 'ifany')
ukb_processed$BMI = ukb_recoded$`03_BMI.0.0`

## CHD coronary heart disease 
# table(ukb_recoded$`04_Angina_age.0.0`, useNA = 'ifany')
# table(ukb_recoded$`04_HeartAttack_age.0.0`, useNA = 'ifany')

## Smoking
table(ukb_recoded$`02_Smoking_status.0.0`, useNA = 'ifany')
ukb_processed$Smoking = as.character(ukb_recoded$`02_Smoking_status.0.0`)

ukb_processed$Smoking = replace(ukb_processed$Smoking, 
                                !ukb_processed$Smoking %in% c('Previous', 'Current'),
                                'Other')

table(ukb_processed$Smoking, useNA = 'ifany')

#################################################################################
# Merge outcome
## Positive Control - Asthma
case_asthma <- output_final_asthma[,c("eid","case")] |> dplyr::rename(case_asthma=case)
ukb_processed <- left_join(ukb_processed, case_asthma)

## Negative Control - Hot drink temp
ukb_processed <- left_join(ukb_processed, NegativeControl)

## Outcome - Dementia
outcome_AD <- output_final_AD[,c("eid","case")] |> dplyr::rename(case_AD = case)
outcome_VD <- output_final_VD[,c("eid","case")] |> dplyr::rename(case_VD = case)
outcome_OD <- output_final_OD[,c("eid","case")] |> dplyr::rename(case_OD = case)

data_full = ukb_processed |>
  left_join(outcome_AD, join_by(eid)) |>
  left_join(outcome_VD, join_by(eid)) |> 
  left_join(outcome_OD, join_by(eid))

## Outcome - exclusion - other types - records 502366
exclusion_criteria.other_types <- (data_full$case_AD + data_full$case_VD + data_full$case_OD) <= 1
### Visualise other cases in each control group
data_full[exclusion_criteria.other_types==FALSE,c("case_AD","case_VD","case_OD")]
### assign NA values to [control] people who are not diagnosis for target subtypes (i.e. AD) but diagnosis for other subtypes 
data_full$case_AD[exclusion_criteria.other_types==FALSE & data_full$case_AD == 0] = NA
data_full$case_VD[exclusion_criteria.other_types==FALSE & data_full$case_VD == 0] = NA
data_full$case_OD[exclusion_criteria.other_types==FALSE & data_full$case_OD == 0] = NA
### check values changed
data_full[exclusion_criteria.other_types==FALSE,c("case_AD","case_VD","case_OD")]

## Define the all cases
data_full$case_all <- as.numeric(data_full$case_AD==1 | data_full$case_VD==1 | data_full$case_OD==1)

## exclude self-report - records 502334
output_final_self_report_dementia <- output_final_self_report_dementia |> dplyr::rename(case_self_report=case)
data_full = data_full |> left_join(output_final_self_report_dementia, join_by(eid))
exlusion_criteria.self_report <- rowSums(data_full[data_full$case_self_report==1,c("case_AD","case_VD","case_OD")], na.rm = TRUE)==0 ## no people having three NA, so na.rm make senes.
exlusion_criteria.self_report <- names(exlusion_criteria.self_report[exlusion_criteria.self_report])
### exclude self-report people who are control or NA for all three subtypes

data_final <- data_full[!(row.names(data_full) %in% exlusion_criteria.self_report), ]

## exclude white ethnicity - records 472542
data_final <- data_final[data_final$ethnicity_white,]

#################################################################################
# save processed data
saveRDS(data_final, "../Data_Extraction/02_Cleaned_Data/Cleaned_ukb_with_outcome.rds")


