
project_root <- normalizePath(Sys.getenv("ITH_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(gtsummary)
library(writexl)
library(zoo)  
library(dplyr)
library(imputeTS)
library(mice)
library(tidyr)

train <- read.csv(file.path(project_root, "data", "clinical", "clinical_data_training.csv"), fileEncoding = "GBK")
valid <- read_excel(file.path(project_root, "data", "clinical", "clinical_data_validation_imputed.xlsx"))

normalize_source_columns <- function(data) {
  source_names <- c("CR", "ID1", "ID2")
  target_names <- c("CR", "Imaging_ID", "Admission_ID")
  for (i in seq_along(source_names)) {
    if (source_names[i] %in% names(data)) names(data)[names(data) == source_names[i]] <- target_names[i]
  }
  data
}
train <- normalize_source_columns(train)
valid <- normalize_source_columns(valid)


train$Cohort <- "Training"
valid$Cohort <- "Validation"
sum(is.na(valid))


valid_filled <- valid %>%
  mutate(
    NLR = as.numeric(NLR)  
  ) 


valid_filled$NLR <- na_kalman(valid_filled$NLR)



valid_filled <- valid_filled %>%
  mutate(
    ITHGroup = as.factor(ITHGroup),
    Sex = as.factor(Sex),
    Tumor_site = as.factor(Tumor_site),
    Pathological_type = as.factor(Pathological_type),
    Histologic.grade_level = as.factor(Histologic.grade_level),
    CEA_level = as.factor(CEA_level),
    CA19_9_level = as.factor(CA19_9_level),
    Nstage = as.factor(Nstage),
    Regimen = as.factor(Regimen),
    Tstage = as.factor(Tstage),
    CR = as.factor(CR),
    Cohort = as.factor(Cohort)
  )
valid_filled$Age <- as.numeric(valid_filled$Age)
valid_filled$NLR <- as.numeric(valid_filled$NLR)
valid_filled$CEA <- as.numeric(valid_filled$CEA)
valid_filled$CA199 <- as.numeric(valid_filled$CA199)
valid_filled$ANC <- as.numeric(valid_filled$ANC)
valid_filled$ALC <- as.numeric(valid_filled$ALC)

valid_filled_for_imputation <- valid_filled %>%
  select(CEA_level, CA19_9_level)  


valid_filled_mice <- mice(valid_filled_for_imputation, method = "polyreg", m = 5, seed = 500)


valid_filled_imputed <- complete(valid_filled_mice, 1)  


valid_filled$CEA_level <- valid_filled_imputed$CEA_level
valid_filled$CA19_9_level <- valid_filled_imputed$CA19_9_level


summary(valid_filled)  

train <- train %>%
  mutate(
    ITHGroup = as.factor(ITHGroup),
    Sex = as.factor(Sex),
    Tumor_site = as.factor(Tumor_site),
    Pathological_type = as.factor(Pathological_type),
    Histologic.grade_level = as.factor(Histologic.grade_level),
    CEA_level = as.factor(CEA_level),
    CA19_9_level = as.factor(CA19_9_level),
    Nstage = as.factor(Nstage),
    Regimen = as.factor(Regimen),
    Tstage = as.factor(Tstage),
    CR = as.factor(CR),
    Cohort = as.factor(Cohort)
  )
train$Age <- as.numeric(train$Age)
train$NLR <- as.numeric(train$NLR)
train$CEA <- as.numeric(train$CEA)
train$CA199 <- as.numeric(train$CA199)
train$Imaging_ID <- as.character(train$Imaging_ID)
train$Admission_ID <- as.character(train$Admission_ID)

valid_filled <- valid_filled %>%
  mutate(
    ITHGroup = as.factor(ITHGroup),
    Sex = as.factor(Sex),
    Tumor_site = as.factor(Tumor_site),
    Pathological_type = as.factor(Pathological_type),
    Histologic.grade_level = as.factor(Histologic.grade_level),
    CEA_level = as.factor(CEA_level),
    CA19_9_level = as.factor(CA19_9_level),
    Nstage = as.factor(Nstage),
    Regimen = as.factor(Regimen),
    Tstage = as.factor(Tstage),
    CR = as.factor(CR),
    Cohort = as.factor(Cohort)
  )
valid_filled$Age <- as.numeric(valid_filled$Age)
valid_filled$NLR <- as.numeric(valid_filled$NLR)
valid_filled$CEA <- as.numeric(valid_filled$CEA)


data_all <- bind_rows(train, valid_filled)


library(gtsummary)
library(dplyr)


table1 <- data_all %>%
  select(
    Cohort,
    ITHGroup, Age, Sex, Tumor_site, Pathological_type,
    Histologic.grade_level, CEA_level, CA19_9_level,
    Nstage, NLR, Regimen, Tstage, CR
  ) %>%
  tbl_summary(
    by = Cohort,          
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",  
      all_categorical() ~ "{n} ({p}%)"               
    ),
    digits = all_continuous() ~ 2  
  ) %>%
  add_overall() %>%  
  add_p(
    test = list(
      all_continuous() ~ "wilcox.test",  
      all_categorical() ~ "chisq.test"   
    ),
    
    pvalue_fun = function(x) style_pvalue(x, digits = 3)  
  ) %>%
  modify_header(label = "**Variable**") %>%  
  bold_labels()  

table1_df <- as_tibble(table1)
write_xlsx(table1_df, file.path(results_dir, "Table1_Baselinenew.xlsx"))


write.csv(valid_filled,file.path(results_dir, "clinical_validation_completed.csv"),fileEncoding = "GB18030")



table_tumor_site <- table(data_all$Tumor_site, data_all$Cohort)


fisher_test_result <- fisher.test(table_tumor_site)


fisher_test_result$p.value




regimen_table <- matrix(
  c(47, 148,19, 60),  
  nrow = 2,             
  ncol = 2,             
  dimnames = list(
    Regimen = c("high", "normal"),
    Cohort = c("Train", "valid")
  )
)


chisq_result <- chisq.test(regimen_table, correct = TRUE)


print(regimen_table)       
print(chisq_result)        
