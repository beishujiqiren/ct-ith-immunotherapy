
project_root <- normalizePath(Sys.getenv("ITH_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(car)
library(broom)
library(openxlsx)

file_path <- file.path(project_root, "data", "clinical", "clinical_data_training.csv")

raw <- read.csv(
  file_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "GB18030"
)


names(raw)

dat <- raw[, c(
  2,   
  3,   
  7,   
  8,   
  10  

  
)]

names(dat) <- c(
  "CR",
  "ITH",
  "Tumor_site",
  "Pathological_type",
  "Histologic_grade"
  

)

dat <- dat %>%
  mutate(
    across(
      everything(),
      ~ if (is.character(.x)) trimws(.x) else .x
    )
  )

cat("\n===== Variable levels =====\n")

for (v in names(dat)) {
  cat("\n", v, ":\n")
  print(table(dat[[v]], useNA = "ifany"))
}

dat$CR <- ifelse(dat$CR == "CR", 1, 0)

table(dat$CR)


dat$ITH <- factor(dat$ITH)
dat$ITH <- relevel(dat$ITH, ref = "Low")



dat$Tumor_site <- factor(dat$Tumor_site)
dat$Tumor_site <- relevel(dat$Tumor_site, ref = "Left-sided")




dat$Pathological_type <- factor(dat$Pathological_type)

non_mucinous_level <- levels(dat$Pathological_type)[
  !grepl(
    "Mucin",
    levels(dat$Pathological_type),
    ignore.case = TRUE
  )
]

dat$Pathological_type <- relevel(
  dat$Pathological_type,
  ref = non_mucinous_level[1]
)




dat$Histologic_grade <- factor(dat$Histologic_grade)

grade_ref <- levels(dat$Histologic_grade)[
  grepl(
    "Moderate|well",
    levels(dat$Histologic_grade),
    ignore.case = TRUE
  )
]

dat$Histologic_grade <- relevel(
  dat$Histologic_grade,
  ref = grade_ref[1]
)





model_formula <- CR ~
  ITH +
  Tumor_site +
  Pathological_type +
  Histologic_grade 
  
 

model_data <- dat %>%
  filter(
    complete.cases(
      ITH,
      Tumor_site,
      Pathological_type,
      Histologic_grade,
      
      
      
      
      CR
    )
  )


cat("\nFinal model N =", nrow(model_data), "\n")
cat("CR =", sum(model_data$CR == 1), "\n")
cat("Non-CR =", sum(model_data$CR == 0), "\n")


fit_main <- glm(
  model_formula,
  data = model_data,
  family = binomial()
)

summary(fit_main)



main_result <- broom::tidy(
  fit_main,
  exponentiate = TRUE,
  conf.int = TRUE
)

main_result <- main_result %>%
  mutate(
    OR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      estimate,
      conf.low,
      conf.high
    )
  )

print(main_result)



vif_result <- car::vif(fit_main)

print(vif_result)



if (is.matrix(vif_result)) {
  
  Table_SX <- data.frame(
    Variable = rownames(vif_result),
    GVIF = vif_result[, "GVIF"],
    Df = vif_result[, "Df"],
    Adjusted_GVIF = vif_result[, ncol(vif_result)],
    row.names = NULL
  )
  
  Table_SX <- Table_SX %>%
    mutate(
      GVIF = round(GVIF, 3),
      Adjusted_GVIF = round(Adjusted_GVIF, 3)
    )
  
  max_vif <- max(
    Table_SX$Adjusted_GVIF,
    na.rm = TRUE
  )
  
} else {
  
  Table_SX <- data.frame(
    Variable = names(vif_result),
    VIF = as.numeric(vif_result)
  )
  
  Table_SX$VIF <- round(
    Table_SX$VIF,
    3
  )
  
  max_vif <- max(
    Table_SX$VIF,
    na.rm = TRUE
  )
}


cat(
  "\nMaximum adjusted GVIF/VIF =",
  round(max_vif, 2),
  "\n"
)

print(Table_SX)











mucinous_flag <- grepl(
  "Mucin",
  model_data$Pathological_type,
  ignore.case = TRUE
)

cat(
  "\nMucinous adenocarcinoma:",
  sum(mucinous_flag),
  "\n"
)






sensitivity_data <- model_data[
  !mucinous_flag,
]


sensitivity_data <- droplevels(
  sensitivity_data
)

cat(
  "Sensitivity analysis N =",
  nrow(sensitivity_data),
  "\n"
)

cat(
  "CR =",
  sum(sensitivity_data$CR == 1),
  "\n"
)

cat(
  "Non-CR =",
  sum(sensitivity_data$CR == 0),
  "\n"
)









sensitivity_formula <- CR ~
  ITH +
  Tumor_site +
  Histologic_grade +
  CEA +
  
  Tstage +
  Nstage


fit_sensitivity <- glm(
  sensitivity_formula,
  data = sensitivity_data,
  family = binomial()
)

summary(fit_sensitivity)






Table_SY <- broom::tidy(
  fit_sensitivity,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    Variable = term,
    
    `OR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      estimate,
      conf.low,
      conf.high
    ),
    
    `P value` = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  )


print(Table_SY)






ITH_sensitivity <- broom::tidy(
  fit_sensitivity,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    grepl(
      "ITH.*High",
      term,
      ignore.case = TRUE
    )
  )

print(ITH_sensitivity)






if (nrow(ITH_sensitivity) == 1) {
  
  cat(
    "\nSensitivity analysis result:\n"
  )
  
  cat(
    "OR =",
    round(ITH_sensitivity$estimate, 2),
    ", 95% CI ",
    round(ITH_sensitivity$conf.low, 2),
    "–",
    round(ITH_sensitivity$conf.high, 2),
    ", P =",
    ifelse(
      ITH_sensitivity$p.value < 0.001,
      "<0.001",
      round(ITH_sensitivity$p.value, 3)
    ),
    "\n",
    sep = ""
  )
}






output_file <- file.path(
  results_dir,
  "Supplementary_Table_SX_SY.xlsx"
)


wb <- createWorkbook()



addWorksheet(
  wb,
  "Table SX - VIF"
)

writeData(
  wb,
  "Table SX - VIF",
  Table_SX
)



addWorksheet(
  wb,
  "Table SY - No mucinous"
)

writeData(
  wb,
  "Table SY - No mucinous",
  Table_SY
)



addWorksheet(
  wb,
  "Main model check"
)

writeData(
  wb,
  "Main model check",
  main_result
)


saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)


cat(
  "\nSaved to:\n",
  output_file,
  "\n"
)









library(dplyr)
library(broom)
library(flextable)
library(officer)







dat2 <- dat %>%
  filter(
    !grepl(
      "Mucin|黏液",
      Pathological_type,
      ignore.case = TRUE
    )
  )





dat2$CR <- ifelse(dat2$CR == "1", 1, 0)

table(dat2$CR)
dat2$ITH <- factor(
  dat2$ITH,
  levels = c("Low", "High")
)

dat2$Tumor_site <- factor(dat2$Tumor_site)
dat2$Tumor_site <- relevel(
  dat2$Tumor_site,
  ref = "Left-sided"
)





dat2$Histologic_grade <- factor(dat2$Histologic_grade)

grade_ref <- levels(dat2$Histologic_grade)[
  grepl(
    "Moderate|well",
    levels(dat2$Histologic_grade),
    ignore.case = TRUE
  )
]

dat2$Histologic_grade <- relevel(
  dat2$Histologic_grade,
  ref = grade_ref[1]
)


dat2$CEA <- factor(dat2$CEA)

cea_ref <- levels(dat2$CEA)[
  grepl(
    "normal|≤5|<=5",
    levels(dat2$CEA),
    ignore.case = TRUE
  )
]

dat2$CEA <- relevel(
  dat2$CEA,
  ref = cea_ref[1]
)












dat2$Tstage <- factor(dat2$Tstage)

t_ref <- levels(dat2$Tstage)[
  grepl(
    "T2/3|T2.*3",
    levels(dat2$Tstage),
    ignore.case = TRUE
  )
]

dat2$Tstage <- relevel(
  dat2$Tstage,
  ref = t_ref[1]
)


dat2$Nstage <- factor(dat2$Nstage)

if ("N+" %in% levels(dat2$Nstage)) {
  dat2$Nstage <- relevel(
    dat2$Nstage,
    ref = "N+"
  )
}






analysis_data <- dat2 %>%
  filter(
    complete.cases(
      CR,
      ITH,
      Tumor_site,
      Histologic_grade,
      CEA,
      
      Tstage,
      Nstage
    )
  )

cat("Sensitivity analysis N =", nrow(analysis_data), "\n")






fit <- glm(
  CR ~
    ITH +
    Tumor_site +
    Histologic_grade +
    CEA +
    
    Tstage +
    Nstage,
  data = analysis_data,
  family = binomial()
)

result <- broom::tidy(
  fit,
  exponentiate = TRUE,
  conf.int = TRUE
)






get_result <- function(pattern) {
  
  x <- result %>%
    filter(
      grepl(
        pattern,
        term,
        ignore.case = TRUE
      )
    )
  
  if (nrow(x) == 0) {
    return(c("", ""))
  }
  
  OR_CI <- sprintf(
    "%.2f (%.2f–%.2f)",
    x$estimate[1],
    x$conf.low[1],
    x$conf.high[1]
  )
  
  P <- ifelse(
    x$p.value[1] < 0.001,
    "<0.001",
    sprintf("%.3f", x$p.value[1])
  )
  
  c(OR_CI, P)
}






get_n <- function(variable, level) {
  
  sum(
    analysis_data[[variable]] == level,
    na.rm = TRUE
  )
}



grade_levels <- levels(analysis_data$Histologic_grade)
cea_levels   <- levels(analysis_data$CEA)
site_levels  <- levels(droplevels(analysis_data$Tumor_site))
t_levels     <- levels(analysis_data$Tstage)
n_levels     <- levels(analysis_data$Nstage)






tab <- data.frame(
  Variables = c(
    
    "ITH",
    "    Low",
    "    High",
    
    "Tumor site",
    paste0("    ", site_levels[1]),
    paste0("    ", site_levels[2]),
    
    "Histologic grade",
    paste0("    ", grade_levels[1]),
    paste0("    ", grade_levels[2]),
    
    "CEA",
    paste0("    ", cea_levels[1]),
    paste0("    ", cea_levels[2]),
    
    "T stage",
    paste0("    ", t_levels[1]),
    paste0("    ", t_levels[2]),
    
    "N stage",
    paste0("    ", n_levels[1]),
    paste0("    ", n_levels[2])
  ),
  
  N = c(
    
    "",
    get_n("ITH", "Low"),
    get_n("ITH", "High"),
    
    "",
    get_n("Tumor_site", site_levels[1]),
    get_n("Tumor_site", site_levels[2]),
    
    "",
    get_n("Histologic_grade", grade_levels[1]),
    get_n("Histologic_grade", grade_levels[2]),
    
    "",
    get_n("CEA", cea_levels[1]),
    get_n("CEA", cea_levels[2]),
    
    "",
    get_n("Tstage", t_levels[1]),
    get_n("Tstage", t_levels[2]),
    
    "",
    get_n("Nstage", n_levels[1]),
    get_n("Nstage", n_levels[2])
  ),
  
  `OR (95% CI)` = c(
    
    "",
    "Reference",
    get_result("^ITHHigh")[1],
    
    "",
    "Reference",
    get_result("^Tumor_site")[1],
    
    "",
    "Reference",
    get_result("Histologic_grade")[1],
    
    "",
    "Reference",
    get_result("^CEA")[1],
    
    "",
    "Reference",
    get_result("^Tstage")[1],
    
    "",
    "Reference",
    get_result("^Nstage")[1]
  ),
  
  `P value` = c(
    
    "",
    "",
    get_result("^ITHHigh")[2],
    
    "",
    "",
    get_result("^Tumor_site")[2],
    
    "",
    "",
    get_result("Histologic_grade")[2],
    
    "",
    "",
    get_result("^CEA")[2],
    
    "",
    "",
    get_result("^Tstage")[2],
    
    "",
    "",
    get_result("^Nstage")[2]
  ),
  
  check.names = FALSE
)


print(tab)






ft <- flextable(tab)


ft <- border_remove(ft)


line_thick <- fp_border(
  color = "black",
  width = 1.2
)

line_thin <- fp_border(
  color = "black",
  width = 0.6
)


ft <- hline_top(
  ft,
  border = line_thick,
  part = "header"
)


ft <- hline_bottom(
  ft,
  border = line_thin,
  part = "header"
)


ft <- hline_bottom(
  ft,
  border = line_thick,
  part = "body"
)


ft <- font(
  ft,
  fontname = "Times New Roman",
  part = "all"
)

ft <- fontsize(
  ft,
  size = 10,
  part = "all"
)


ft <- bold(
  ft,
  part = "header"
)


ft <- align(
  ft,
  j = 1,
  align = "left",
  part = "all"
)

ft <- align(
  ft,
  j = 2:4,
  align = "center",
  part = "all"
)


group_rows <- c(
  1, 4, 7, 10, 13, 16
)

ft <- bold(
  ft,
  i = group_rows,
  j = 1,
  bold = TRUE
)


p_numeric <- suppressWarnings(
  as.numeric(tab$`P value`)
)

sig_rows <- which(
  !is.na(p_numeric) &
    p_numeric < 0.05
)

ft <- bold(
  ft,
  i = sig_rows,
  j = 4,
  bold = TRUE
)


ft <- autofit(ft)






doc <- read_docx()

doc <- body_add_par(
  doc,
  "Supplementary Table SY. Sensitivity analysis of factors associated with complete response after excluding patients with mucinous adenocarcinoma",
  style = "Normal"
)

doc <- body_add_flextable(
  doc,
  ft
)

doc <- body_add_par(
  doc,
  "Note: ORs were estimated using multivariable logistic regression. Left-sided tumor location was used as the reference category. ITH, intratumoral heterogeneity; CEA, carcinoembryonic antigen.",
  style = "Normal"
)

output_path <- file.path(
  results_dir,
  "Supplementary_Table_SY.docx"
)

print(
  doc,
  target = output_path
)

cat(
  "\nTable saved to:\n",
  output_path,
  "\n"
)

