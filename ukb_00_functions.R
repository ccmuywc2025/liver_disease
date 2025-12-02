# This script contains all custom functions required for UKB external validation analysis
# Function: Define functions for data processing, statistical analysis, table generation, and plotting
# Input: Data frames and parameters
# Output: Processed data, statistical result lists, chart files

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tableone))
suppressPackageStartupMessages(library(pROC))

# 1. Outlier processing - Winsorize method
winsorize_variable <- function(x, probs = c(0.01, 0.99)) {
  if (all(is.na(x))) return(x)
  q <- quantile(x, probs, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

# 2. Standardized variable creation
standardize_variable <- function(x) {
  if (all(is.na(x))) return(x)
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

# 3. Summarize Cox regression results
summarize_cox_results <- function(cox_results_list, model_type = "univariate") {
  results_df <- data.frame(
    Exposure = character(),
    HR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    P_value = numeric(),
    N_events = numeric(),
    N_total = numeric(),
    Model_type = character()
  )
  
  for (var in names(cox_results_list)) {
    summary_obj <- cox_results_list[[var]]
    hr <- summary_obj$conf.int[1, "exp(coef)"]
    ci_lower <- summary_obj$conf.int[1, "lower .95"]
    ci_upper <- summary_obj$conf.int[1, "upper .95"]
    p_value <- summary_obj$coef[1, 5]
    n_events <- summary_obj$nevent
    model_label <- ifelse(tolower(model_type) == "multivariate", "Multivariate", "Univariate")
    
    results_df <- rbind(results_df, data.frame(
      Exposure = var,
      HR = hr,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      P_value = p_value,
      N_events = n_events,
      N_total = summary_obj$n,
      Model_type = model_label
    ))
  }
  
  return(results_df)
}

# 4. Summarize disease-specific results
summarize_disease_specific_results <- function(disease_results_list) {
  all_results <- data.frame()
  
  for (disease_name in names(disease_results_list)) {
    disease_results <- disease_results_list[[disease_name]]
    
    # Univariate results
    if (!is.null(disease_results$univariate)) {
      univ_df <- summarize_cox_results(disease_results$univariate, "univariate")
      univ_df$Disease_type <- disease_name
      univ_df$Model_type <- "Univariate"
      all_results <- rbind(all_results, univ_df)
    }
    
    # Multivariate results
    if (!is.null(disease_results$multivariate)) {
      multiv_df <- summarize_cox_results(disease_results$multivariate, "multivariate")
      multiv_df$Disease_type <- disease_name
      multiv_df$Model_type <- "Multivariate"
      all_results <- rbind(all_results, multiv_df)
    }
  }
  
  return(all_results)
}

# 5. Summarize Logistic regression results
summarize_logistic_results <- function(res_list) {
  out <- data.frame(
    Exposure = character(),
    OR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    P_value = numeric(),
    N_total = numeric(),
    Model_type = character()
  )
  for (nm in names(res_list)) {
    sm <- res_list[[nm]]
    or <- exp(sm$coefficients[2, "Estimate"])
    ci <- exp(sm$confint[2, c("2.5 %", "97.5 %")])
    pv <- sm$coefficients[2, "Pr(>|z|)"]
    nt <- sm$N_total
    mt <- sm$Model_type
    out <- rbind(out, data.frame(
      Exposure = nm,
      OR = or,
      CI_lower = ci[1],
      CI_upper = ci[2],
      P_value = pv,
      N_total = nt,
      Model_type = mt
    ))
  }
  out
}

# 6. Execute Logistic regression analysis set
fit_logistic_set <- function(df, outcome, exposures, covariates, model_label) {
  univ <- list(); multiv <- list()
  for (var in exposures) {
    if (var %in% names(df)) {
      f_univ <- as.formula(paste(outcome, "~", var))
      m_univ <- glm(f_univ, data = df, family = binomial())
      s_univ <- summary(m_univ)
      ci_univ <- suppressWarnings(confint(m_univ))
      nt_univ <- sum(complete.cases(df[, c(outcome, var)]))
      s_univ$confint <- ci_univ
      s_univ$N_total <- nt_univ
      s_univ$Model_type <- paste(model_label, "Univariate")
      univ[[var]] <- s_univ
      cov_terms <- paste(covariates[covariates %in% names(df)], collapse = " + ")
      if (nzchar(cov_terms)) {
        f_mult <- as.formula(paste(outcome, "~", var, "+", cov_terms))
      } else {
        f_mult <- as.formula(paste(outcome, "~", var))
      }
      m_mult <- glm(f_mult, data = df, family = binomial())
      s_mult <- summary(m_mult)
      ci_mult <- suppressWarnings(confint(m_mult))
      nt_mult <- sum(complete.cases(df[, c(outcome, var, covariates)]))
      s_mult$confint <- ci_mult
      s_mult$N_total <- nt_mult
      s_mult$Model_type <- paste(model_label, "Multivariate")
      multiv[[var]] <- s_mult
    }
  }
  list(
    univariate = summarize_logistic_results(univ),
    multivariate = summarize_logistic_results(multiv)
  )
}

# 7. Create UKB baseline table
create_ukb_baseline_table <- function(data) {
  candidate_vars <- c(
    "p21022",
    "p31",
    "p21001_i0",
    "p20116_i0_3cat",
    "p1558_i0_3cat",
    "p6138_i0_3cat",
    "p20118_i0_urban_rural",
    "cti_ukb",
    "egdr_ukb",
    "tyg_ukb",
    "tyg_whtr_ukb"
  )
  vars <- candidate_vars[candidate_vars %in% names(data)]
  factor_vars <- c("p31", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural")
  factor_vars <- factor_vars[factor_vars %in% names(data)]
  strata_var <- "incident_cld"
  tab <- CreateTableOne(vars = vars, strata = strata_var, data = data, factorVars = factor_vars, test = TRUE)
  out <- print(tab, showAllLevels = TRUE, formatOptions = list(big.mark = ","), printToggle = FALSE)
  dir.create("results/ukb", showWarnings = FALSE, recursive = TRUE)
  write.csv(out, "results/ukb/baseline_table.csv")
  wb_baseline <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_baseline, "Baseline_Table")
  openxlsx::writeData(wb_baseline, "Baseline_Table", out)
  openxlsx::addStyle(wb_baseline, "Baseline_Table", openxlsx::createStyle(halign = "center"), rows = 1, cols = 1:ncol(out))
  openxlsx::setColWidths(wb_baseline, "Baseline_Table", cols = 1:ncol(out), widths = "auto")
  openxlsx::saveWorkbook(wb_baseline, "results/ukb/baseline_table.xlsx", overwrite = TRUE)
  list(table = tab, output = out)
}

# 8. Create MRI baseline table
create_mri_baseline_table <- function(data, mri_var, year_label) {
  use_data <- data[!is.na(data[[mri_var]]), ]
  candidate_vars <- c(
    "p21022",
    "p31",
    "p21001_i0",
    "p20116_i0_3cat",
    "p1558_i0_3cat",
    "p6138_i0_3cat",
    "p20118_i0_urban_rural",
    "cti_ukb",
    "egdr_ukb",
    "tyg_ukb",
    "tyg_whtr_ukb"
  )
  vars <- candidate_vars[candidate_vars %in% names(use_data)]
  factor_vars <- c("p31", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural")
  factor_vars <- factor_vars[factor_vars %in% names(use_data)]
  tab <- CreateTableOne(vars = vars, strata = mri_var, data = use_data, factorVars = factor_vars, test = TRUE)
  out <- print(tab, showAllLevels = TRUE, formatOptions = list(big.mark = ","), printToggle = FALSE)
  dir.create("results/ukb", showWarnings = FALSE, recursive = TRUE)
  csv_path <- paste0("results/ukb/mri_baseline_table_", year_label, ".csv")
  xlsx_path <- paste0("results/ukb/mri_baseline_table_", year_label, ".xlsx")
  write.csv(out, csv_path)
  wb_mri_base <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_mri_base, paste0("MRI", year_label, "Baseline_Table"))
  openxlsx::writeData(wb_mri_base, paste0("MRI", year_label, "Baseline_Table"), out)
  openxlsx::addStyle(wb_mri_base, paste0("MRI", year_label, "Baseline_Table"), openxlsx::createStyle(halign = "center"), rows = 1, cols = 1:ncol(out))
  openxlsx::setColWidths(wb_mri_base, paste0("MRI", year_label, "Baseline_Table"), cols = 1:ncol(out), widths = "auto")
  openxlsx::saveWorkbook(wb_mri_base, xlsx_path, overwrite = TRUE)
  list(table = tab, output = out)
}

# 9. Compute MRI ROC and AUC
compute_mri_roc_auc <- function(df, outcome, exposures, covariates, year_label) {
  res <- data.frame(
    Year = character(),
    Exposure = character(),
    Model = character(),
    AUC = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    N_total = integer(),
    stringsAsFactors = FALSE
  )
  dir.create("results/ukb", showWarnings = FALSE, recursive = TRUE)
  label_map <- c(
    cti_ukb = "CTI Index",
    egdr_ukb = "eGDR",
    tyg_ukb = "TyG Index",
    tyg_whtr_ukb = "TyG-WHtR"
  )
  for (var in exposures) {
    if (!(var %in% names(df))) next
    df_univ <- df[complete.cases(df[, c(outcome, var)]), ]
    if (nrow(df_univ) > 0) {
      m_univ <- glm(as.formula(paste(outcome, "~", var)), data = df_univ, family = binomial())
      phat_univ <- predict(m_univ, type = "response")
      roc_univ <- pROC::roc(df_univ[[outcome]], phat_univ, direction = "auto", quiet = TRUE)
      ci_univ <- as.numeric(pROC::ci.auc(roc_univ))
      tiff_path <- paste0("results/ukb/roc_", year_label, "_", var, "_univ.tiff")
      tiff(tiff_path, width = 2400, height = 2000, res = 300)
      plot(roc_univ, print.auc = TRUE)
      dev.off()
      res <- rbind(res, data.frame(
        Year = year_label,
        Exposure = ifelse(var %in% names(label_map), unname(label_map[var]), var),
        Model = "Univariate",
        AUC = as.numeric(pROC::auc(roc_univ)),
        CI_lower = ci_univ[1],
        CI_upper = ci_univ[3],
        N_total = nrow(df_univ),
        stringsAsFactors = FALSE
      ))
    }
    cov_terms <- covariates[covariates %in% names(df)]
    df_mult <- df[complete.cases(df[, c(outcome, var, cov_terms)]), ]
    if (nrow(df_mult) > 0) {
      if (length(cov_terms) > 0) {
        f_str <- paste(outcome, "~", var, "+", paste(cov_terms, collapse = " + "))
      } else {
        f_str <- paste(outcome, "~", var)
      }
      m_mult <- glm(as.formula(f_str), data = df_mult, family = binomial())
      phat_mult <- predict(m_mult, type = "response")
      roc_mult <- pROC::roc(df_mult[[outcome]], phat_mult, direction = "auto", quiet = TRUE)
      ci_mult <- as.numeric(pROC::ci.auc(roc_mult))
      tiff_path <- paste0("results/ukb/roc_", year_label, "_", var, "_multiv.tiff")
      tiff(tiff_path, width = 2400, height = 2000, res = 300)
      plot(roc_mult, print.auc = TRUE)
      dev.off()
      res <- rbind(res, data.frame(
        Year = year_label,
        Exposure = ifelse(var %in% names(label_map), unname(label_map[var]), var),
        Model = "Multivariate",
        AUC = as.numeric(pROC::auc(roc_mult)),
        CI_lower = ci_mult[1],
        CI_upper = ci_mult[3],
        N_total = nrow(df_mult),
        stringsAsFactors = FALSE
      ))
    }
  }
  csv_path <- paste0("results/ukb/mri_roc_auc_results_", year_label, ".csv")
  write.csv(res, csv_path, row.names = FALSE, fileEncoding = "gbk")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, paste0("MRI", year_label))
  openxlsx::writeData(wb, paste0("MRI", year_label), res)
  openxlsx::addStyle(wb, paste0("MRI", year_label), openxlsx::createStyle(halign = "center"), rows = 1, cols = 1:ncol(res))
  openxlsx::setColWidths(wb, paste0("MRI", year_label), cols = 1:ncol(res), widths = "auto")
  openxlsx::saveWorkbook(wb, paste0("results/ukb/mri_roc_auc_results_", year_label, ".xlsx"), overwrite = TRUE)
  res
}
