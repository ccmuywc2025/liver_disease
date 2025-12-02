# UKB External Validation Cox Regression Analysis Script
# Function: Perform Cox regression analysis (overall liver disease risk and specific liver disease types)
# Input: Processed RDS data file
# Output: Cox regression results (Excel and RDS)


source("R/code_release/ukb/ukb_00_functions.R")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(openxlsx))

# Read Data
df <- readRDS("data/processed/ukb_validation_dataset.rds")

# Define Analysis Variables
exposure_vars <- c("cti_ukb_std", "egdr_ukb_std", "tyg_ukb_std", "tyg_whtr_ukb_std")
covariates <- c("p21022", "p31", "p21001_i0", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural")
# Ensure Variables Exist
exposure_vars <- exposure_vars[exposure_vars %in% names(df)]
covariates <- covariates[covariates %in% names(df)]

cat(sprintf("\n=== Cox Regression Analysis Start ===\n"))
cat(sprintf("Sample Size: %d\n", nrow(df)))

# 1) Overall Liver Disease Risk
surv_obj <- Surv(time = df$follow_up_time_years, event = df$incident_cld)

# Univariate Analysis
cox_univariate_results <- list()
for (var in exposure_vars) {
  if (var %in% names(df)) {
    formula <- as.formula(paste("surv_obj ~", var))
    cox_model <- coxph(formula, data = df)
    cox_univariate_results[[var]] <- summary(cox_model)
  }
}

# Multivariate Analysis
cox_multivariate_results <- list()
for (var in exposure_vars) {
  if (var %in% names(df)) {
    covar_terms <- paste(covariates, collapse = " + ")
    formula <- as.formula(paste("surv_obj ~", var, "+", covar_terms))
    cox_model <- coxph(formula, data = df)
    cox_multivariate_results[[var]] <- summary(cox_model)
  }
}

univ_summary <- summarize_cox_results(cox_univariate_results, "univariate")
multiv_summary <- summarize_cox_results(cox_multivariate_results, "multivariate")

if (nrow(univ_summary) > 0) {
  univ_summary$P_adj_BH <- p.adjust(univ_summary$P_value, method = "BH")
}
if (nrow(multiv_summary) > 0) {
  multiv_summary$P_adj_BH <- p.adjust(multiv_summary$P_value, method = "BH")
}

# 2) Specific Liver Disease Type Analysis
liver_disease_types <- list(
  "alcoholic_liver" = "p131658_after",          # K70: Alcoholic liver disease
  "toxic_liver" = "p131660_after",              # K71: Toxic liver disease
  "hepatitis" = "p131662_after",                  # K72: Hepatitis
  "cirrhosis" = "p131664_after",                  # K73: Cirrhosis
  "fibrosis" = "p131666_after",                   # K74: Fibrosis
  "other_liver" = "p131670_after"                   # K76: Other liver diseases
)

disease_specific_results <- list()
for (disease_name in names(liver_disease_types)) {
  event_var <- liver_disease_types[[disease_name]]
  
  if (event_var %in% names(df)) {
    time_var <- paste0(event_var, "_end_time_years")
    event_var_event <- paste0(event_var, "_event")
    
    if (time_var %in% names(df) && event_var_event %in% names(df)) {
      df[[time_var]] <- as.numeric(df[[time_var]])
      surv_obj_specific <- Surv(time = df[[time_var]], event = df[[event_var_event]])
      
      disease_results <- list()
      # Univariate
      for (var in exposure_vars) {
        if (var %in% names(df)) {
          formula <- as.formula(paste("surv_obj_specific ~", var))
          cox_model <- coxph(formula, data = df)
          disease_results$univariate[[var]] <- summary(cox_model)
        }
      }
      # Multivariate
      for (var in exposure_vars) {
        if (var %in% names(df)) {
          covar_terms <- paste(covariates, collapse = " + ")
          formula <- as.formula(paste("surv_obj_specific ~", var, "+", covar_terms))
          cox_model <- coxph(formula, data = df)
          disease_results$multivariate[[var]] <- summary(cox_model)
        }
      }
      disease_specific_results[[disease_name]] <- disease_results
    }
  }
}

disease_summary <- summarize_disease_specific_results(disease_specific_results)
if (nrow(disease_summary) > 0) {
  disease_summary$P_adj_BH <- NA_real_
  grp <- split(seq_len(nrow(disease_summary)), paste(disease_summary$Disease_type, disease_summary$Model_type))
  for (idxs in grp) {
    pv <- disease_summary$P_value[idxs]
    disease_summary$P_adj_BH[idxs] <- p.adjust(pv, method = "BH")
  }
}

# 3) Output Results
output_results_path <- "results/ukb/cox_regression_results.rds"
output_excel_path <- "results/ukb/cox_regression_results.xlsx"
dir.create(dirname(output_results_path), showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  univariate_results = cox_univariate_results,
  multivariate_results = cox_multivariate_results,
  disease_specific_results = disease_specific_results,
  univariate_summary = univ_summary,
  multivariate_summary = multiv_summary,
  disease_specific_summary = disease_summary
), file = output_results_path)

wb <- createWorkbook()

if (nrow(univ_summary) > 0) {
  d <- univ_summary[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Model_Type")
  addWorksheet(wb, "Overall_Univariate")
  writeData(wb, "Overall_Univariate", d)
}
if (nrow(multiv_summary) > 0) {
  d <- multiv_summary[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Model_Type")
  addWorksheet(wb, "Overall_Multivariate")
  writeData(wb, "Overall_Multivariate", d)
}
if (nrow(disease_summary) > 0) {
  d <- disease_summary[, c("Exposure","HR","CI_lower","CI_upper","P_value","P_adj_BH","N_events","N_total","Disease_type","Model_type")]
  names(d) <- c("Exposure", "HR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Events", "N_Total", "Disease_Type", "Model_Type")
  addWorksheet(wb, "Disease_Specific")
  writeData(wb, "Disease_Specific", d)
}

saveWorkbook(wb, output_excel_path, overwrite = TRUE)
cat(sprintf("Cox regression results saved to %s\n", output_excel_path))
