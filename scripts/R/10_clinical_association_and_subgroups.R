



project_root <- normalizePath(Sys.getenv("ITH_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "clinical_association")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(logistf)
})
if (!requireNamespace("car", quietly = TRUE)) {
  stop("Package 'car' is required for GVIF. Install it with install.packages('car').")
}


INPUT_FILE <- file.path(project_root, "data", "clinical", "clinical_data_training.csv")
OUTCOME_COLUMN <- "CR"                 
ITH_COLUMN <- "ITHGroup"                    
CANDIDATE_PREDICTORS <- c(
  "ITHGroup", "Tumor_site_type", "Pathological_type", "Histologic.grade_level"
)
SUBGROUP_VARIABLES <- c("Tumor_site_type", "Pathological_type", "Histologic.grade_level")


UNIVARIABLE_P_THRESHOLD <- 0.10

REFERENCE_PREFERENCES <- list(
  ITHGroup = "Low",
  Sex = "female",
  Tumor_site_type = "Left-sided",
  Pathological_type = "Adenocarcinoma",
  Histologic.grade_level = "moderate"
)

if (!file.exists(INPUT_FILE)) stop("Training clinical file not found: ", INPUT_FILE)
raw <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "GB18030")
if (!OUTCOME_COLUMN %in% names(raw) && "CR" %in% names(raw)) {
  names(raw)[names(raw) == "CR"] <- OUTCOME_COLUMN
}
if (!OUTCOME_COLUMN %in% names(raw)) stop("Outcome column not found: ", OUTCOME_COLUMN)
if (!ITH_COLUMN %in% names(raw)) stop("ITH column not found: ", ITH_COLUMN)

available_predictors <- intersect(CANDIDATE_PREDICTORS, names(raw))
missing_predictors <- setdiff(CANDIDATE_PREDICTORS, names(raw))
if (length(missing_predictors)) {
  stop("Required Methods variables missing from the clinical file: ",
       paste(missing_predictors, collapse = ", "),
       ". Check the names above against the source CSV and edit the configuration.")
}


cr_text <- trimws(as.character(raw[[OUTCOME_COLUMN]]))
cr_lower <- tolower(cr_text)
outcome <- ifelse(cr_lower %in% c("cr", "1", "yes", "true"), 1L,
                  ifelse(cr_lower %in% c("non-cr", "noncr", "0", "no", "false"), 0L, NA_integer_))
if (any(is.na(outcome) & !is.na(cr_text) & nzchar(cr_text))) {
  bad <- unique(cr_text[is.na(outcome) & !is.na(cr_text) & nzchar(cr_text)])
  stop("Unrecognized outcome values: ", paste(bad, collapse = ", "),
       ". Recode these explicitly before analysis.")
}
dat <- raw[, unique(c(available_predictors)), drop = FALSE]
dat$Outcome_CR <- outcome


categorical_predictors <- intersect(
  c("ITHGroup", "Tumor_site_type", "Pathological_type",
    "Histologic.grade_level"), names(dat)
)
for (v in categorical_predictors) {
  dat[[v]] <- factor(trimws(as.character(dat[[v]])))
  pref <- REFERENCE_PREFERENCES[[v]]
  if (!is.null(pref) && length(levels(dat[[v]]))) {
    lev <- levels(dat[[v]])
    hit <- which(tolower(lev) == tolower(pref))
    if (!length(hit) && identical(v, "Pathological_type"))
      hit <- which(!grepl("mucin", lev, ignore.case = TRUE))
    if (!length(hit) && identical(v, "Histologic.grade_level"))
      hit <- which(grepl("moderate|well", lev, ignore.case = TRUE))
    if (length(hit)) {
      dat[[v]] <- relevel(dat[[v]], ref = lev[hit[1]])
    } else {
      message("Reference preference '", pref, "' not found for ", v,
              "; using first factor level: ", lev[1])
    }
  }
}

if (!all(c("Low", "High") %in% levels(dat[[ITH_COLUMN]]))) {
  stop("ITHGroup must contain exactly the expected labels 'Low' and 'High'. Found: ",
       paste(levels(dat[[ITH_COLUMN]]), collapse = ", "))
}
dat[[ITH_COLUMN]] <- relevel(dat[[ITH_COLUMN]], ref = "Low")


usable <- available_predictors[vapply(available_predictors, function(v) {
  x <- dat[[v]]
  length(unique(x[!is.na(x)])) > 1
}, logical(1))]
if (!ITH_COLUMN %in% usable) stop("ITHGroup has no usable variation.")
write.csv(data.frame(candidate = available_predictors,
                     included = available_predictors %in% usable),
          file.path(results_dir, "predictor_columns_used.csv"), row.names = FALSE)


get_lrt_p <- function(fit, variable) {
  tab <- drop1(fit, test = "LRT")
  if (!variable %in% rownames(tab)) return(NA_real_)
  pcol <- grep("Pr\\(>Chi\\)", colnames(tab), value = TRUE)
  if (!length(pcol)) return(NA_real_)
  as.numeric(tab[variable, pcol[1]])
}

uni_summary <- list()
uni_terms <- list()
for (v in usable) {
  d <- dat[complete.cases(dat[, c("Outcome_CR", v), drop = FALSE]), c("Outcome_CR", v), drop = FALSE]
  if (nrow(d) < 5 || length(unique(d$Outcome_CR)) < 2 || length(unique(d[[v]])) < 2) next
  fit <- tryCatch(glm(reformulate(v, response = "Outcome_CR"), data = d,
                      family = binomial()), error = function(e) NULL)
  if (is.null(fit)) next
  p_global <- get_lrt_p(fit, v)
  uni_summary[[v]] <- data.frame(
    variable = v, n = nrow(d), events_CR = sum(d$Outcome_CR == 1),
    global_LR_P = p_global, stringsAsFactors = FALSE
  )
  cf <- summary(fit)$coefficients
  if (nrow(cf) > 1) {
    ci <- suppressMessages(confint.default(fit))
    term_names <- rownames(cf)
    term_names <- term_names[term_names != "(Intercept)"]
    if (length(term_names)) {
      uni_terms[[v]] <- data.frame(
        variable = v, term = term_names,
        OR = exp(coef(fit)[term_names]),
        CI_low = exp(ci[term_names, 1]), CI_high = exp(ci[term_names, 2]),
        P = cf[term_names, "Pr(>|z|)"], n = nrow(d),
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
}
uni_summary_df <- if (length(uni_summary)) do.call(rbind, uni_summary) else data.frame()
uni_terms_df <- if (length(uni_terms)) do.call(rbind, uni_terms) else data.frame()
write.csv(uni_summary_df, file.path(results_dir, "univariable_logistic_global_tests.csv"), row.names = FALSE)
write.csv(uni_terms_df, file.path(results_dir, "univariable_logistic_coefficients.csv"), row.names = FALSE)


selected <- if (nrow(uni_summary_df)) {
  uni_summary_df$variable[!is.na(uni_summary_df$global_LR_P) &
                            uni_summary_df$global_LR_P < UNIVARIABLE_P_THRESHOLD]
} else character()
if (!length(selected)) stop("No predictor met the prespecified univariable P < 0.10 threshold.")
write.csv(data.frame(selected_predictor = selected),
          file.path(results_dir, "multivariable_predictors_selected.csv"), row.names = FALSE)


main_vars <- selected
main_dat <- dat[complete.cases(dat[, c("Outcome_CR", main_vars), drop = FALSE]),
                c("Outcome_CR", main_vars), drop = FALSE]
main_fit <- glm(reformulate(main_vars, response = "Outcome_CR"), data = main_dat,
                family = binomial())
main_cf <- summary(main_fit)$coefficients
main_ci <- suppressMessages(confint.default(main_fit))
main_terms <- rownames(main_cf)
main_terms <- main_terms[main_terms != "(Intercept)"]
main_table <- data.frame(
  term = main_terms, OR = exp(coef(main_fit)[main_terms]),
  CI_low = exp(main_ci[main_terms, 1]), CI_high = exp(main_ci[main_terms, 2]),
  P = main_cf[main_terms, "Pr(>|z|)"], n = nrow(main_dat),
  events_CR = sum(main_dat$Outcome_CR == 1), stringsAsFactors = FALSE,
  row.names = NULL
)
write.csv(main_table, file.path(results_dir, "multivariable_logistic.csv"), row.names = FALSE)

vif_result <- tryCatch(car::vif(main_fit), error = function(e) NULL)
if (!is.null(vif_result)) {
  if (is.matrix(vif_result)) {
    vif_df <- data.frame(variable = rownames(vif_result), vif_result,
                         row.names = NULL, check.names = FALSE)
    if (all(c("GVIF", "Df") %in% names(vif_df)))
      vif_df$GVIF_adjusted <- vif_df$GVIF^(1 / (2 * vif_df$Df))
  } else {
    vif_df <- data.frame(variable = names(vif_result), GVIF = as.numeric(vif_result),
                         row.names = NULL)
  }
  write.csv(vif_df, file.path(results_dir, "multivariable_gvif.csv"), row.names = FALSE)
}




ith_numeric <- ifelse(dat[[ITH_COLUMN]] == "High", 1, 0)
spear_vars <- setdiff(intersect(CANDIDATE_PREDICTORS, names(dat)), ITH_COLUMN)
spear_out <- list()
for (v in spear_vars) {
  x <- dat[[v]]
  if (is.factor(x) || is.character(x)) {
    lev <- levels(factor(x))
    if (length(lev) > 2) {
      message("Spearman skipped for ", v, ": >2 unordered levels; define an ordinal coding if justified.")
      next
    }
    x <- as.integer(factor(x, levels = lev))
  } else {
    x <- suppressWarnings(as.numeric(x))
  }
  ok <- is.finite(x) & !is.na(ith_numeric)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2) next
  ct <- suppressWarnings(cor.test(ith_numeric[ok], x[ok], method = "spearman", exact = FALSE))
  spear_out[[v]] <- data.frame(variable = v, n = sum(ok), rho = unname(ct$estimate),
                               P = ct$p.value, stringsAsFactors = FALSE)
}
spear_df <- if (length(spear_out)) do.call(rbind, spear_out) else data.frame()
write.csv(spear_df, file.path(results_dir, "spearman_ith_clinical_correlations.csv"), row.names = FALSE)


subgroups <- intersect(SUBGROUP_VARIABLES, names(dat))
subgroups <- subgroups[subgroups %in% usable]
cc_vars <- unique(c("Outcome_CR", ITH_COLUMN, main_vars, subgroups))
cc_dat <- dat[complete.cases(dat[, cc_vars, drop = FALSE]), cc_vars, drop = FALSE]
subgroup_results <- list()
interaction_results <- list()

firth_ith_row <- function(fit, subgroup, level, n, covariates) {
  coef_names <- names(fit$coefficients)
  effect_term <- coef_names[grepl(paste0("^", ITH_COLUMN), coef_names) &
                              grepl("High", coef_names, ignore.case = TRUE)]
  if (!length(effect_term)) {
    return(data.frame(subgroup = subgroup, level = level, n = n,
                      adjusted_for = paste(covariates, collapse = ";"),
                      OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                      note = "ITH High-vs-Low coefficient not estimable", stringsAsFactors = FALSE))
  }
  term <- effect_term[1]
  data.frame(subgroup = subgroup, level = level, n = n,
             adjusted_for = paste(covariates, collapse = ";"),
             OR = exp(fit$coefficients[term]), CI_low = exp(fit$ci.lower[term]),
             CI_high = exp(fit$ci.upper[term]), P = fit$prob[term],
             note = "", stringsAsFactors = FALSE)
}

for (sg in subgroups) {
  for (lev in levels(droplevels(cc_dat[[sg]]))) {
    d <- droplevels(cc_dat[cc_dat[[sg]] == lev, , drop = FALSE])
    if (length(unique(d$Outcome_CR)) < 2 || length(unique(d[[ITH_COLUMN]])) < 2) {
      subgroup_results[[length(subgroup_results) + 1L]] <- data.frame(
        subgroup = sg, level = lev, n = nrow(d), adjusted_for = NA_character_,
        OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
        note = "Insufficient outcome or ITH-group variation", stringsAsFactors = FALSE)
      next
    }
    covars <- setdiff(main_vars, sg)
    covars <- covars[vapply(covars, function(v) {
      length(unique(d[[v]][!is.na(d[[v]])])) > 1
    }, logical(1))]
    covars <- setdiff(covars, ITH_COLUMN)
    f <- reformulate(c(ITH_COLUMN, covars), response = "Outcome_CR")
    fit <- tryCatch(logistf::logistf(f, data = d, pl = TRUE), error = function(e) e)
    if (inherits(fit, "error")) {
      subgroup_results[[length(subgroup_results) + 1L]] <- data.frame(
        subgroup = sg, level = lev, n = nrow(d), adjusted_for = paste(covars, collapse = ";"),
        OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
        note = conditionMessage(fit), stringsAsFactors = FALSE)
    } else {
      subgroup_results[[length(subgroup_results) + 1L]] <- firth_ith_row(fit, sg, lev, nrow(d), covars)
    }
  }

  
  base_vars <- unique(c(main_vars, sg, ITH_COLUMN))
  base_vars <- base_vars[vapply(base_vars, function(v) {
    v == ITH_COLUMN || length(unique(cc_dat[[v]][!is.na(cc_dat[[v]])])) > 1
  }, logical(1))]
  base_vars <- unique(c(ITH_COLUMN, setdiff(base_vars, ITH_COLUMN)))
  main_formula <- reformulate(base_vars, response = "Outcome_CR")
  interaction_formula <- as.formula(paste(
    "Outcome_CR ~", paste(base_vars, collapse = " + "), "+", ITH_COLUMN, ":", sg
  ))
  fit_int <- tryCatch(logistf::logistf(interaction_formula, data = cc_dat, pl = TRUE),
                      error = function(e) e)
  if (inherits(fit_int, "error")) {
    interaction_results[[sg]] <- data.frame(
      subgroup = sg, n = nrow(cc_dat), LR_chisq = NA_real_, df = NA_real_,
      P_interaction = NA_real_, note = conditionMessage(fit_int), stringsAsFactors = FALSE)
  } else {
    test <- tryCatch(anova(fit_int, formula = as.formula(paste("~", ITH_COLUMN, ":", sg)),
                           method = "nested"), error = function(e) e)
    if (inherits(test, "error")) {
      interaction_results[[sg]] <- data.frame(
        subgroup = sg, n = nrow(cc_dat), LR_chisq = NA_real_, df = NA_real_,
        P_interaction = NA_real_, note = conditionMessage(test), stringsAsFactors = FALSE)
    } else {
      interaction_results[[sg]] <- data.frame(
        subgroup = sg, n = nrow(cc_dat), LR_chisq = as.numeric(test$chisq),
        df = as.numeric(test$df), P_interaction = as.numeric(test$pval),
        note = "", stringsAsFactors = FALSE)
    }
  }
}
subgroup_df <- if (length(subgroup_results)) do.call(rbind, subgroup_results) else data.frame()
interaction_df <- if (length(interaction_results)) do.call(rbind, interaction_results) else data.frame()
write.csv(subgroup_df, file.path(results_dir, "firth_subgroup_estimates.csv"), row.names = FALSE)
write.csv(interaction_df, file.path(results_dir, "firth_interaction_tests.csv"), row.names = FALSE)

cat("Analysis complete. Results saved under: ", results_dir, "\n", sep = "")
cat("Main model predictors selected by univariable LRT P < 0.10: ",
    paste(selected, collapse = ", "), "\n", sep = "")
cat("Subgroup analyses used a common complete-case population (n = ", nrow(cc_dat), ").\n", sep = "")
