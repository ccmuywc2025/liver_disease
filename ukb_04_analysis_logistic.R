# UKB External Validation Logistic Regression and MRI Analysis Script
# Function: Perform Logistic regression for MRI-detected fatty liver (MASLD), generate baseline table, and ROC analysis
# Input: Processed RDS data file
# Output: Logistic regression results, baseline table, ROC curves and results


source("R/code_release/ukb/ukb_00_functions.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tableone))
suppressPackageStartupMessages(library(pROC))

df <- readRDS("data/processed/ukb_validation_dataset.rds")

# Derive MRI Outcome Variables
df <- df %>%
  mutate(
    masld_mri_2019 = if_else(is.na(p40061_i3), as.integer(NA), if_else(p40061_i3 > 5, 1L, 0L)),
    masld_mri_2014 = if_else(is.na(p40061_i2), as.integer(NA), if_else(p40061_i2 > 5, 1L, 0L))
  )

# Define Analysis Variables
exposure_vars_logit <- c("cti_ukb_std", "egdr_ukb_std", "tyg_ukb_std", "tyg_whtr_ukb_std")
covariates_logit <- c("p21022", "p31", "p21001_i0", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural")
exposure_vars_logit <- exposure_vars_logit[exposure_vars_logit %in% names(df)]
covariates_logit <- covariates_logit[covariates_logit %in% names(df)]

# 1. Logistic Regression Analysis
mri2014_results <- fit_logistic_set(df, "masld_mri_2014", exposure_vars_logit, covariates_logit, "MRI2014")
mri2019_results <- fit_logistic_set(df, "masld_mri_2019", exposure_vars_logit, covariates_logit, "MRI2019")

# Adjust P-values
adjust_p_values <- function(res) {
  if (nrow(res) > 0) {
    res$P_adj_BH <- p.adjust(res$P_value, method = "BH")
  }
  return(res)
}

if (!is.null(mri2014_results$univariate)) mri2014_results$univariate <- adjust_p_values(mri2014_results$univariate)
if (!is.null(mri2014_results$multivariate)) mri2014_results$multivariate <- adjust_p_values(mri2014_results$multivariate)
if (!is.null(mri2019_results$univariate)) mri2019_results$univariate <- adjust_p_values(mri2019_results$univariate)
if (!is.null(mri2019_results$multivariate)) mri2019_results$multivariate <- adjust_p_values(mri2019_results$multivariate)

# Save Logistic Results
output_results_path_mri_logit <- "results/ukb/masld_mri_logistic_results.rds"
output_excel_path_mri_logit <- "results/ukb/masld_mri_logistic_results.xlsx"
dir.create(dirname(output_results_path_mri_logit), showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  mri2014_univ = mri2014_results$univariate,
  mri2014_mult = mri2014_results$multivariate,
  mri2019_univ = mri2019_results$univariate,
  mri2019_mult = mri2019_results$multivariate
), file = output_results_path_mri_logit)

wb_mri_logit <- createWorkbook()
save_logit_sheet <- function(wb, data, sheet_name) {
  if (!is.null(data) && nrow(data) > 0) {
    d <- data[, c("Exposure","OR","CI_lower","CI_upper","P_value","P_adj_BH","N_total","Model_type")]
    names(d) <- c("Exposure", "OR", "95%CI_Lower", "95%CI_Upper", "P_Value", "P_adj_BH", "N_Total", "Model_Type")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, d)
  }
}

save_logit_sheet(wb_mri_logit, mri2014_results$univariate, "MRI2014_Univariate")
save_logit_sheet(wb_mri_logit, mri2014_results$multivariate, "MRI2014_Multivariate")
save_logit_sheet(wb_mri_logit, mri2019_results$univariate, "MRI2019_Univariate")
save_logit_sheet(wb_mri_logit, mri2019_results$multivariate, "MRI2019_Multivariate")
saveWorkbook(wb_mri_logit, output_excel_path_mri_logit, overwrite = TRUE)
cat(sprintf("Logistic regression results saved to %s\n", output_excel_path_mri_logit))

# 2. Baseline Table Generation
baseline_table_res <- create_ukb_baseline_table(df)
mri_table_2014 <- create_mri_baseline_table(df, "masld_mri_2014", "2014")
mri_table_2019 <- create_mri_baseline_table(df, "masld_mri_2019", "2019")

# 3. ROC Analysis (Logistic)
roc_auc_2014 <- compute_mri_roc_auc(df, "masld_mri_2014", exposure_vars_logit, covariates_logit, "2014")
roc_auc_2019 <- compute_mri_roc_auc(df, "masld_mri_2019", exposure_vars_logit, covariates_logit, "2019")

cat("MRI Logistic regression and related analysis completed\n")
