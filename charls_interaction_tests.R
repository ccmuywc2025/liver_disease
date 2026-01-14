
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

# Load data
data <- readRDS("data/processed/merged_individual_blood_2020.rds")
data <- data[inw1==1]

# Baseline exclusions
baseline_participants <- data[!is.na(r1livere)]
longitudinal_data <- baseline_participants[r1livere == 0]

# Covariate exclusions
longitudinal_data <- longitudinal_data[!is.na(nl_Blood_weight)]
longitudinal_data <- longitudinal_data[!is.na(r1agey)]
longitudinal_data <- longitudinal_data[!is.na(r1mbmi)]
longitudinal_data <- longitudinal_data[!is.na(r1smoken)]
longitudinal_data <- longitudinal_data[!is.na(r1drinkev)]

# Define outcome
longitudinal_data[, liver_incident := ifelse((r2livere == 1 & !is.na(r2livere)) | 
                                             (r3livere == 1 & !is.na(r3livere)) | 
                                             (r4livere == 1 & !is.na(r4livere)) | 
                                             (r5livere == 1 & !is.na(r5livere)), 1, 0)]

# Define biomarkers and standardize
main_biomarkers <- c("cti_2011", "egdr_2011", "tyg_2011", "tyg_whtr_2011")
biomarker_names <- c("cti_2011"="CTI(2011)", "egdr_2011"="eGDR(2011)", 
                     "tyg_2011"="TyG(2011)", "tyg_whtr_2011"="TyG-WHtR(2011)")

for (bm in main_biomarkers) {
  std <- paste0(bm, "_std")
  v <- longitudinal_data[[bm]]
  idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident)
  mu <- mean(v[idx], na.rm=TRUE)
  sdv <- sd(v[idx], na.rm=TRUE)
  if (is.finite(sdv) && sdv != 0) {
    longitudinal_data[, (std) := (get(bm)-mu)/sdv]
  } else {
    longitudinal_data[, (std) := NA_real_]
  }
}

# Define subgroups variables
longitudinal_data[, age_group := ifelse(r1agey >= 60, 1, 0)]
longitudinal_data[, bmi_group := ifelse(r1mbmi >= 24, 1, 0)]

subgroup_vars <- list(
  "Gender" = "ragender",
  "Age Group" = "age_group",
  "BMI Group" = "bmi_group",
  "Residence" = "h1rural"
)

interaction_results <- data.frame()

# Base covariates
base_covariates_list <- c("r1agey", "ragender", "r1mbmi", "raeducl", "h1rural", 
                          "marital", "r1smoken", "r1drinkev", "r1chronic_group")

for (biomarker in main_biomarkers) {
  std_var <- paste0(biomarker, "_std")
  
  for (subgroup_name in names(subgroup_vars)) {
    subgroup_var <- subgroup_vars[[subgroup_name]]
    
    # Filter data
    analysis_data <- longitudinal_data[!is.na(get(std_var)) & !is.na(liver_incident) & !is.na(get(subgroup_var))]
    
    # Adjust covariates: remove the subgroup variable itself from adjustment
    covariates <- base_covariates_list
    if (subgroup_var == "ragender") covariates <- setdiff(covariates, "ragender")
    if (subgroup_var == "age_group") covariates <- setdiff(covariates, "r1agey")
    if (subgroup_var == "bmi_group") covariates <- setdiff(covariates, "r1mbmi")
    if (subgroup_var == "h1rural") covariates <- setdiff(covariates, "h1rural")
    
    # Model without interaction
    formula_no_int <- as.formula(paste("liver_incident ~", std_var, "+", subgroup_var, "+", 
                                       paste(covariates, collapse = " + ")))
    model_no_int <- glm(formula_no_int, data = analysis_data, family = binomial())
    
    # Model with interaction
    formula_int <- as.formula(paste("liver_incident ~", std_var, "*", subgroup_var, "+", 
                                    paste(covariates, collapse = " + ")))
    model_int <- glm(formula_int, data = analysis_data, family = binomial())
    
    # LRT Test
    lrt <- anova(model_no_int, model_int, test = "Chisq")
    p_interaction <- lrt$`Pr(>Chi)`[2]
    
    interaction_results <- rbind(interaction_results, data.frame(
      biomarker = biomarker,
      biomarker_name = biomarker_names[biomarker],
      subgroup = subgroup_name,
      interaction_p_lrt = p_interaction,
      stringsAsFactors = FALSE
    ))
  }
}

# Save results
dir.create("results/charls", showWarnings = FALSE, recursive = TRUE)
write.csv(interaction_results, "results/charls/interaction_test_results.csv", row.names = FALSE)
cat("Interaction tests completed. Results saved to results/charls/interaction_test_results.csv\n")
