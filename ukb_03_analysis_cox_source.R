# UKB External Validation Cox Regression Analysis Script (Source Confirmed Outcomes)
# Function: Filter liver disease events by source (Death register/Hospital admissions) and perform Cox regression
# Input: Processed RDS data file
# Output: Source-confirmed Cox regression results (Excel and RDS)


source("R/code_release/ukb/ukb_00_functions.R")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(dplyr))

df <- readRDS("data/processed/ukb_validation_dataset.rds")
data_lock_date <- max(df$follow_up_end_date, na.rm = TRUE) # Approximate lock date or recalculate

# Recalculate data_lock_date (if precision needed)
# Infer from df, or use max(df$p53_i0 + df$follow_up_time_days, na.rm=TRUE)
# In original script, data_lock_date was a global variable. Here we recalculate.
# Use max of follow_up_end_date for consistency
if (!exists("data_lock_date")) {
  data_lock_date <- max(as.Date(df$follow_up_end_date), na.rm = TRUE)
}

# Define Confirmation Sources
allowed_sources <- tolower(c(
  "Death register only",
  "Death register and other source(s)",
  "Hospital admissions data only",
  "Hospital admissions data and other source(s)"
))

pairs <- list(
  c("p131658","p131659"),
  c("p131660","p131661"),
  c("p131662","p131663"),
  c("p131664","p131665"),
  c("p131666","p131667"),
  c("p131668","p131669"),
  c("p131670","p131671")
)

for (pr in pairs) {
  ev <- pr[1]; src <- pr[2]
  if (ev %in% names(df) && src %in% names(df)) {
    nm <- paste0(ev, "_after_src")
    df[[nm]] <- dplyr::if_else(
      !is.na(df[[ev]]) & !is.na(df[["p53_i0"]]) & (df[[ev]] >= df[["p53_i0"]]) &
        (tolower(trimws(df[[src]])) %in% allowed_sources),
      df[[ev]],
      as.Date(NA)
    )
  }
}

after_src_cols <- intersect(paste0(c("p131658","p131660","p131662","p131664","p131666","p131668","p131670"), "_after_src"), names(df))
df <- df %>%
  mutate(
    cld_event_date_src = {
      if (length(after_src_cols) == 0) {
        as.Date(rep(NA, n()))
      } else {
        cols <- across(all_of(after_src_cols))
        out <- suppressWarnings(do.call(pmin, c(cols, na.rm = TRUE)))
        if (!inherits(out, "Date")) {
          out <- as.Date(out, origin = "1970-01-01")
        }
        na_all <- rowSums(!is.na(cols)) == 0L
        out[na_all] <- as.Date(NA)
        out
      }
    }
  ) %>%
  mutate(
    incident_cld_src = if_else(is.na(cld_event_date_src), 0L, 1L)
  ) %>%
  mutate(
    follow_up_end_date_src = if_else(incident_cld_src == 1L & !is.na(cld_event_date_src), cld_event_date_src, data_lock_date),
    follow_up_time_days_src = if_else(is.na(p53_i0) | is.na(follow_up_end_date_src), as.numeric(NA), pmax(0, as.numeric(difftime(follow_up_end_date_src, p53_i0, units = "days")))),
    follow_up_time_years_src = if_else(is.na(follow_up_time_days_src), as.numeric(NA), follow_up_time_days_src / 365.25)
  ) %>%
  mutate(
    across(
      all_of(after_src_cols),
      ~ if_else(is.na(.x), 0L, 1L),
      .names = "{.col}_event"
    )
  ) %>%
  mutate(
    across(
      all_of(after_src_cols),
      ~ if_else(!is.na(.x), .x, data_lock_date),
      .names = "{.col}_end"
    )
  ) %>%
  mutate(
    across(
      all_of(paste0(after_src_cols, "_end")),
      ~ if_else(is.na(p53_i0) | is.na(.x), as.numeric(NA), pmax(0, as.numeric(difftime(.x, p53_i0, units = "days")))),
      .names = "{.col}_time_days"
    )
  ) %>%
  mutate(
    across(
      all_of(paste0(after_src_cols, "_end")),
      ~ as.numeric(if_else(is.na(p53_i0) | is.na(.x), NA_real_, pmax(0, as.numeric(difftime(.x, p53_i0, units = "days"))) / 365.25)),
      .names = "{.col}_time_years"
    )
  )

# Define Analysis Variables
exposure_vars_src <- c("cti_ukb_std","egdr_ukb_std", "tyg_ukb_std", "tyg_whtr_ukb_std")
covariates_src <- c("p21022", "p31", "p21001_i0", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural")
exposure_vars_src <- exposure_vars_src[exposure_vars_src %in% names(df)]
covariates_src <- covariates_src[covariates_src %in% names(df)]

# 1) Overall Liver Disease Risk (Source Confirmed)
surv_obj_src <- Surv(time = df$follow_up_time_years_src, event = df$incident_cld_src)

# Univariate
cox_univariate_results_src <- list()
for (var in exposure_vars_src) {
  if (var %in% names(df)) {
    formula <- as.formula(paste("surv_obj_src ~", var))
    cox_model <- coxph(formula, data = df)
    cox_univariate_results_src[[var]] <- summary(cox_model)
  }
}

# Multivariate
cox_multivariate_results_src <- list()
for (var in exposure_vars_src) {
  if (var %in% names(df)) {
    covar_terms <- paste(covariates_src, collapse = " + ")
    formula <- as.formula(paste("surv_obj_src ~", var, "+", covar_terms))
    cox_model <- coxph(formula, data = df)
    cox_multivariate_results_src[[var]] <- summary(cox_model)
  }
}

univ_summary_src <- summarize_cox_results(cox_univariate_results_src, "univariate")
multiv_summary_src <- summarize_cox_results(cox_multivariate_results_src, "multivariate")

if (nrow(univ_summary_src) > 0) {
  univ_summary_src$P_adj_BH <- p.adjust(univ_summary_src$P_value, method = "BH")
}
if (nrow(multiv_summary_src) > 0) {
  multiv_summary_src$P_adj_BH <- p.adjust(multiv_summary_src$P_value, method = "BH")
}

# 2) Specific Liver Disease Type (Source Confirmed)
liver_disease_types_src <- list(
  "alcoholic_liver" = "p131658_after_src",
  "toxic_liver" = "p131660_after_src",
  "hepatitis" = "p131662_after_src",
  "cirrhosis" = "p131664_after_src",
  "fibrosis" = "p131666_after_src",
  "other_liver" = "p131670_after_src"
)

disease_specific_results_src <- list()
for (disease_name in names(liver_disease_types_src)) {
  event_var <- liver_disease_types_src[[disease_name]]
  if (event_var %in% names(df)) {
    time_var <- paste0(event_var, "_end_time_years")
    event_var_event <- paste0(event_var, "_event")
    if (time_var %in% names(df) && event_var_event %in% names(df)) {
      df[[time_var]] <- as.numeric(df[[time_var]])
      surv_obj_specific <- Surv(time = df[[time_var]], event = df[[event_var_event]])
      
      disease_results <- list()
      # Univariate
      for (var in exposure_vars_src) {
        if (var %in% names(df)) {
          formula <- as.formula(paste("surv_obj_specific ~", var))
          cox_model <- coxph(formula, data = df)
          disease_results$univariate[[var]] <- summary(cox_model)
        }
      }
      # Multivariate
      for (var in exposure_vars_src) {
        if (var %in% names(df)) {
          covar_terms <- paste(covariates_src, collapse = " + ")
          formula <- as.formula(paste("surv_obj_specific ~", var, "+", covar_terms))
          cox_model <- coxph(formula, data = df)
          disease_results$multivariate[[var]] <- summary(cox_model)
        }
      }
      disease_specific_results_src[[disease_name]] <- disease_results
    }
  }
}

disease_summary_src <- summarize_disease_specific_results(disease_specific_results_src)
if (nrow(disease_summary_src) > 0) {
  disease_summary_src$P_adj_BH <- NA_real_
  grp_src <- split(seq_len(nrow(disease_summary_src)), paste(disease_summary_src$Disease_type, disease_summary_src$Model_type))
  for (idxs in grp_src) {
    pv <- disease_summary_src$P_value[idxs]
    disease_summary_src$P_adj_BH[idxs] <- p.adjust(pv, method = "BH")
  }
}

# 3) Output Results
output_results_path_src <- "results/ukb/cox_regression_results_src.rds"
output_excel_path_src <- "results/ukb/cox_regression_results_src.xlsx"
dir.create(dirname(output_results_path_src), showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  univariate_results = cox_univariate_results_src,
  multivariate_results = cox_multivariate_results_src,
  disease_specific_results = disease_specific_results_src,
  univariate_summary = univ_summary_src,
  multivariate_summary = multiv_summary_src,
  disease_specific_summary = disease_summary_src
), file = output_results_path_src)

wb_src <- createWorkbook()
if (nrow(univ_summary_src) > 0) {
  d <- univ_summary_src[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Model_Type")
  addWorksheet(wb_src, "Overall_Univariate_src")
  writeData(wb_src, "Overall_Univariate_src", d)
}
if (nrow(multiv_summary_src) > 0) {
  d <- multiv_summary_src[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Model_Type")
  addWorksheet(wb_src, "Overall_Multivariate_src")
  writeData(wb_src, "Overall_Multivariate_src", d)
}
if (nrow(disease_summary_src) > 0) {
  d <- disease_summary_src[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Disease_type","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Disease_Type", "Model_Type")
  addWorksheet(wb_src, "Disease_Specific_src")
  writeData(wb_src, "Disease_Specific_src", d)
}
saveWorkbook(wb_src, output_excel_path_src, overwrite = TRUE)
cat(sprintf("Source-confirmed Cox regression results saved to %s\n", output_excel_path_src))
