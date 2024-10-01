# Load required library
library(dplyr)
library(parallel)
library(table1)
library(ggplot2)

library(gtsummary)

#no_cores=detectCores()-1
#cl <- makeCluster(no_cores) 
#################################################################################
# Load Data
#setwd("/rds/general/project/hda_23-24/live/TDS/group1/Scripts")

## load processed ukb dataset with lables
data_full <- readRDS("../Data_Extraction/02_Cleaned_Data/Cleaned_ukb_with_PRS_new.rds")
#################################################################################
# Build Table 1

data_full$case_all <- 
  factor(data_full$case_all, 
         levels=c(1,0),
         labels=c("Case", # Reference
                  "Control"))
table(data_full$case_all, useNA = 'ifany')

data_full$case_AD <- 
  factor(data_full$case_AD, 
         levels=c(1,0),
         labels=c("Case", # Reference
                  "Control"))
table(data_full$case_AD, useNA = 'ifany')

data_full$case_VD <- 
  factor(data_full$case_VD, 
         levels=c(1,0),
         labels=c("Case", # Reference
                  "Control"))
table(data_full$case_VD, useNA = 'ifany')

data_full$case_OD <- 
  factor(data_full$case_OD, 
         levels=c(1,0),
         labels=c("Case", # Reference
                  "Control"))
table(data_full$case_OD, useNA = 'ifany')

table_all <- data_full %>% 
  select(case_all, age_now, gender, edu_level, #ethnicity_white, 
         SBP, TotalChol, BMI, Smoking, HotDrinkTemp_Binary) %>%
  tbl_summary(
    by = case_all,
    statistic = list(all_categorical() ~ "{n}    ({p}%)",
                     age_now     ~ "{mean} ({sd})",
                     SBP         ~ "{mean} ({sd})",
                     TotalChol   ~ "{mean} ({sd})",
                     BMI         ~ "{mean} ({sd})"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(gender               ~ "categorical",
                edu_level            ~ "categorical",
                #ethnicity_white      ~ "categorical",
                Smoking              ~ "categorical",
                HotDrinkTemp_Binary  ~ "categorical",
                age_now              ~ "continuous",
                SBP                  ~ "continuous",
                TotalChol            ~ "continuous",
                BMI                  ~ "continuous"),
    label = list(age_now    ~ "Current Age (Years)",
                 gender  ~ "Gender",
                 edu_level     ~ "Education Level (if higher education)",
                 #ethnicity_white     ~ "Ethnicity (if white)",
                 SBP   ~ "Systolic Blood Presure (mmHg)",
                 TotalChol  ~ "Total Cholesterol (mmol/l)",
                 BMI  ~ "Body Mass Index (kg/m2)",
                 Smoking ~ "Smoking",
                 HotDrinkTemp_Binary ~ "If take hot drinks") 
  ) %>%
  modify_header(
    label = "**Variable**",
    # The following adds the % to the column total label
    # <br> is the location of a line break
    all_stat_cols() ~ "**{level}**<br>N = {n} ({style_percent(p, digits=1)}%)"
  ) %>%
  add_p() %>%
  bold_labels()  %>%
  # Include an "overall" column
  add_overall(
    last = FALSE,
    # The ** make it bold
    col_label = "**All participants**<br>N = {N}"
  )

table_AD <- data_full %>% 
  select(case_AD, age_now, gender, edu_level, #ethnicity_white, 
         SBP, TotalChol, BMI, Smoking, HotDrinkTemp_Binary) %>%
  tbl_summary(
    by = case_AD,
    statistic = list(all_categorical() ~ "{n}    ({p}%)",
                     age_now     ~ "{mean} ({sd})",
                     SBP         ~ "{mean} ({sd})",
                     TotalChol   ~ "{mean} ({sd})",
                     BMI         ~ "{mean} ({sd})"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(gender               ~ "categorical",
                edu_level            ~ "categorical",
                #ethnicity_white      ~ "categorical",
                Smoking              ~ "categorical",
                HotDrinkTemp_Binary  ~ "categorical",
                age_now              ~ "continuous",
                SBP                  ~ "continuous",
                TotalChol            ~ "continuous",
                BMI                  ~ "continuous"),
    label = list(age_now    ~ "Current Age (Years)",
                 gender  ~ "Gender",
                 edu_level     ~ "Education Level (if higher education)",
                 #ethnicity_white     ~ "Ethnicity (if white)",
                 SBP   ~ "Systolic Blood Presure (mmHg)",
                 TotalChol  ~ "Total Cholesterol (mmol/l)",
                 BMI  ~ "Body Mass Index (kg/m2)",
                 Smoking ~ "Smoking",
                 HotDrinkTemp_Binary ~ "If take hot drinks") 
  ) %>%
  modify_header(
    label = "**Variable**",
    # The following adds the % to the column total label
    # <br> is the location of a line break
    all_stat_cols() ~ "**{level}**<br>N = {n} ({style_percent(p, digits=1)}%)"
  ) %>%
  add_p() %>%
  bold_labels()

table_VD <- data_full %>% 
  select(case_VD, age_now, gender, edu_level, #ethnicity_white, 
         SBP, TotalChol, BMI, Smoking, HotDrinkTemp_Binary) %>%
  tbl_summary(
    by = case_VD,
    statistic = list(all_categorical() ~ "{n}    ({p}%)",
                     age_now     ~ "{mean} ({sd})",
                     SBP         ~ "{mean} ({sd})",
                     TotalChol   ~ "{mean} ({sd})",
                     BMI         ~ "{mean} ({sd})"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(gender               ~ "categorical",
                edu_level            ~ "categorical",
                #ethnicity_white      ~ "categorical",
                Smoking              ~ "categorical",
                HotDrinkTemp_Binary  ~ "categorical",
                age_now              ~ "continuous",
                SBP                  ~ "continuous",
                TotalChol            ~ "continuous",
                BMI                  ~ "continuous"),
    label = list(age_now    ~ "Current Age (Years)",
                 gender  ~ "Gender",
                 edu_level     ~ "Education Level (if higher education)",
                 #ethnicity_white     ~ "Ethnicity (if white)",
                 SBP   ~ "Systolic Blood Presure (mmHg)",
                 TotalChol  ~ "Total Cholesterol (mmol/l)",
                 BMI  ~ "Body Mass Index (kg/m2)",
                 Smoking ~ "Smoking",
                 HotDrinkTemp_Binary ~ "If take hot drinks") 
  ) %>%
  modify_header(
    label = "**Variable**",
    # The following adds the % to the column total label
    # <br> is the location of a line break
    all_stat_cols() ~ "**{level}**<br>N = {n} ({style_percent(p, digits=1)}%)"
  ) %>%
  add_p() %>%
  bold_labels()

table_OD <- data_full %>% 
  select(case_OD, age_now, gender, edu_level, #ethnicity_white, 
         SBP, TotalChol, BMI, Smoking, HotDrinkTemp_Binary) %>%
  tbl_summary(
    by = case_OD,
    statistic = list(all_categorical() ~ "{n}    ({p}%)",
                     age_now     ~ "{mean} ({sd})",
                     SBP         ~ "{mean} ({sd})",
                     TotalChol   ~ "{mean} ({sd})",
                     BMI         ~ "{mean} ({sd})"),
    digits = list(all_continuous()  ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(gender               ~ "categorical",
                edu_level            ~ "categorical",
                #ethnicity_white      ~ "categorical",
                Smoking              ~ "categorical",
                HotDrinkTemp_Binary  ~ "categorical",
                age_now              ~ "continuous",
                SBP                  ~ "continuous",
                TotalChol            ~ "continuous",
                BMI                  ~ "continuous"),
    label = list(age_now    ~ "Current Age (Years)",
                 gender  ~ "Gender",
                 edu_level     ~ "Education Level (if higher education)",
                 #ethnicity_white     ~ "Ethnicity (if white)",
                 SBP   ~ "Systolic Blood Presure (mmHg)",
                 TotalChol  ~ "Total Cholesterol (mmol/l)",
                 BMI  ~ "Body Mass Index (kg/m2)",
                 Smoking ~ "Smoking",
                 HotDrinkTemp_Binary ~ "If take hot drinks") 
  ) %>%
  modify_header(
    label = "**Variable**",
    # The following adds the % to the column total label
    # <br> is the location of a line break
    all_stat_cols() ~ "**{level}**<br>N = {n} ({style_percent(p, digits=1)}%)"
  ) %>%
  add_p() %>%
  bold_labels()


# merge tables
tbl_merge <-
  tbl_merge(
    tbls = list(table_all, table_AD, table_VD, table_OD),
    tab_spanner = c("**All types**", "**Alzheimer's Disease**", "**Vascular Dementia**", "**Other Dementia**")
  )

tbl_merge


pdf(file = "../plots/table1.pdf", width =20, height = 15) # The height of the plot in inches
tbl_merge
dev.off()

tbl_merge %>%
  as_gt() %>%
  gt::gtsave(filename = "../plots/table1.pdf")