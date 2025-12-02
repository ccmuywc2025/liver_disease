suppressPackageStartupMessages(library(data.table))
data <- readRDS("data/processed/merged_individual_blood_2020.rds")
data <- data[inw1==1]
baseline_participants <- data[!is.na(r1livere)]
longitudinal_data <- baseline_participants[r1livere == 0]
longitudinal_data <- longitudinal_data[!is.na(nl_Blood_weight)]
longitudinal_data <- longitudinal_data[!is.na(r1agey)]
longitudinal_data <- longitudinal_data[!is.na(r1mbmi)]
longitudinal_data <- longitudinal_data[!is.na(r1smoken)]
longitudinal_data <- longitudinal_data[!is.na(r1drinkev)]
longitudinal_data[, liver_incident := ifelse((r2livere == 1 & !is.na(r2livere)) | (r3livere == 1 & !is.na(r3livere)) | (r4livere == 1 & !is.na(r4livere)) | (r5livere == 1 & !is.na(r5livere)), 1, 0)]
# 1. Interaction analysis (Continuous variables)
# Create a list to store results
interaction_results_cont <- list()

# Define biomarkers to analyze
biomarkers <- c("z_score_log_tg_hdl", "z_score_log_tyg", "z_score_log_tyg_bmi", "z_score_log_tyg_wc")

# Define covariates
covariates_interaction <- c("age_base", "gender", "bmi_base", "smoking_status", "drinking_status", 
                            "hypertension", "diabetes", "dyslipidemia", "cvd_history")

# Loop through each biomarker
for (biom in biomarkers) {
  # Adjust covariates (remove components included in the biomarker)
  covars <- covariates_interaction
  if (grepl("bmi", biom)) covars <- setdiff(covars, "bmi_base")
  if (grepl("dyslipidemia", biom)) covars <- setdiff(covars, "dyslipidemia")
  if (grepl("tg", biom)) covars <- setdiff(covars, "dyslipidemia")
  
  # Interaction terms to test
  interactions <- c("gender", "age_base", "bmi_base", "diabetes")
  
  for (inter_var in interactions) {
    # Skip if interaction variable is part of the biomarker
    if (inter_var == "bmi_base" && grepl("bmi", biom)) next
    if (inter_var == "diabetes" && grepl("tyg", biom)) next # TyG relates to glucose
    
    # Construct formula
    formula_str <- paste("Surv(time_to_event, status) ~", biom, "*", inter_var, "+", 
                         paste(setdiff(covars, inter_var), collapse = " + "))
    
    fit <- coxph(as.formula(formula_str), data = df)
    
    # Extract P-value for interaction
    p_val <- summary(fit)$coefficients[nrow(summary(fit)$coefficients), "Pr(>|z|)"]
    
    interaction_results_cont[[paste(biom, inter_var, sep = "_")]] <- data.frame(
      Biomarker = biom,
      Interaction_Var = inter_var,
      P_Interaction = p_val
    )
  }
}

# Combine results
res_interaction_cont <- do.call(rbind, interaction_results_cont)
write_xlsx(list(Interaction_Continuous = res_interaction_cont), "results/charls_interaction_continuous.xlsx")

# 2. Subgroup analysis (Forest plot data)
# Define subgroups
subgroup_vars <- list(
  "Gender" = list(var = "ragender", levels = c("Male" = 1, "Female" = 2)),
  "Age" = list(var = "age_group", levels = c("<60 years" = 0, ">=60 years" = 1)),
  "BMI Group" = list(var = "bmi_group", levels = c("Normal" = 0, "Overweight/Obese" = 1)),
  "Smoking Status" = list(var = "smoking_status", levels = c("Never" = 0, "Ever/Current" = 1)),
  "Drinking Status" = list(var = "drinking_status", levels = c("Never" = 0, "Ever/Current" = 1)),
  "Hypertension" = list(var = "hypertension", levels = c("No" = 0, "Yes" = 1)),
  "Diabetes" = list(var = "diabetes", levels = c("No" = 0, "Yes" = 1)),
  "Dyslipidemia" = list(var = "dyslipidemia", levels = c("No" = 0, "Yes" = 1))
)

subgroup_results <- list()

for (biom in biomarkers) {
  for (sub_name in names(subgroup_vars)) {
    sub_info <- subgroup_vars[[sub_name]]
    var_name <- sub_info$var
    
    for (level_name in names(sub_info$levels)) {
      level_val <- sub_info$levels[level_name]
      
      # Subset data
      sub_df <- df[df[[var_name]] == level_val, ]
      
      # Adjust covariates
      covars <- covariates_interaction
      # Remove stratification variable from covariates
      covars <- setdiff(covars, var_name)
      # Remove components
      if (grepl("bmi", biom)) covars <- setdiff(covars, "bmi_base")
      if (grepl("dyslipidemia", biom)) covars <- setdiff(covars, "dyslipidemia")
      if (grepl("tg", biom)) covars <- setdiff(covars, "dyslipidemia")
      
      # Cox model
      formula_str <- paste("Surv(time_to_event, status) ~", biom, "+", paste(covars, collapse = " + "))
      
      tryCatch({
        fit <- coxph(as.formula(formula_str), data = sub_df)
        res <- summary(fit)$coefficients[biom, ]
        
        subgroup_results[[paste(biom, sub_name, level_name, sep = "_")]] <- data.frame(
          Biomarker = biom,
          Subgroup = sub_name,
          Level = level_name,
          HR = exp(res["coef"]),
          Lower = exp(res["coef"] - 1.96 * res["se(coef)"]),
          Upper = exp(res["coef"] + 1.96 * res["se(coef)"]),
          P_Value = res["Pr(>|z|)"],
          N = nrow(sub_df)
        )
      }, error = function(e) {
        message(paste("Error in", biom, sub_name, level_name, ":", e$message))
      })
    }
  }
}

res_subgroup <- do.call(rbind, subgroup_results)
write_xlsx(list(Interaction_Stratified = res_subgroup), "results/charls_subgroup_forest.xlsx")

cat("Interaction and subgroup analysis completed.\n")

