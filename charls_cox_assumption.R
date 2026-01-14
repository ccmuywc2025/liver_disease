# Cox Proportional Hazards Assumption Tests
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

# Load data
data <- readRDS("data/processed/merged_individual_blood_2020.rds")
data <- data[inw1==1]
baseline_participants <- data[!is.na(r1livere)]
longitudinal_data <- baseline_participants[r1livere == 0]

# Apply inclusion criteria
longitudinal_data <- longitudinal_data[!is.na(nl_Blood_weight)]
longitudinal_data <- longitudinal_data[!is.na(r1agey)]
longitudinal_data <- longitudinal_data[!is.na(r1mbmi)]
longitudinal_data <- longitudinal_data[!is.na(r1smoken)]
longitudinal_data <- longitudinal_data[!is.na(r1drinkev)]

# Define incident liver disease
longitudinal_data[, liver_incident := ifelse((r2livere == 1 & !is.na(r2livere)) | 
                                             (r3livere == 1 & !is.na(r3livere)) | 
                                             (r4livere == 1 & !is.na(r4livere)) | 
                                             (r5livere == 1 & !is.na(r5livere)), 1, 0)]

# Calculate follow-up time (Mid-point imputation)
longitudinal_data[, follow_up_time := { 
  time <- rep(0.5, .N)
  incident_cases <- which(!is.na(liver_incident) & liver_incident == 1)
  for (i in incident_cases) { 
    if (!is.na(r2livere[i]) && r2livere[i] == 1) { 
      time[i] <- (0 + 2) / 2 
    } else if (!is.na(r3livere[i]) && r3livere[i] == 1) { 
      last_visit_time <- 0
      if (!is.na(r2livere[i]) && r2livere[i] == 0) last_visit_time <- 2
      time[i] <- (last_visit_time + 4) / 2
    } else if (!is.na(r4livere[i]) && r4livere[i] == 1) { 
      last_visit_time <- 0
      if (!is.na(r3livere[i]) && r3livere[i] == 0) last_visit_time <- 4
      else if (!is.na(r2livere[i]) && r2livere[i] == 0) last_visit_time <- 2
      time[i] <- (last_visit_time + 7) / 2
    } else if (!is.na(r5livere[i]) && r5livere[i] == 1) { 
      last_visit_time <- 0
      if (!is.na(r4livere[i]) && r4livere[i] == 0) last_visit_time <- 7
      else if (!is.na(r3livere[i]) && r3livere[i] == 0) last_visit_time <- 4
      else if (!is.na(r2livere[i]) && r2livere[i] == 0) last_visit_time <- 2
      time[i] <- (last_visit_time + 9) / 2
    } 
  }
  non_incident <- which(is.na(liver_incident) | liver_incident == 0)
  for (i in non_incident) { 
    if (!is.na(r5livere[i])) { time[i] <- 9 } 
    else if (!is.na(r4livere[i])) { time[i] <- 7 } 
    else if (!is.na(r3livere[i])) { time[i] <- 4 } 
    else if (!is.na(r2livere[i])) { time[i] <- 2 } 
  }
  time 
}]
longitudinal_data[follow_up_time <= 0, follow_up_time := 0.5]

# Standardize biomarkers
main_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)")

for (bm in main_biomarkers) { 
  std <- paste0(bm,"_std")
  v <- longitudinal_data[[bm]]
  idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident)
  mu <- mean(v[idx], na.rm=TRUE)
  sdv <- sd(v[idx], na.rm=TRUE)
  if (is.finite(sdv) && sdv != 0) longitudinal_data[, (std) := (get(bm)-mu)/sdv] else longitudinal_data[, (std) := NA_real_] 
}

# Perform Cox PH assumption tests
cox_results <- data.frame()

for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(longitudinal_data))) next
  
  # Filter data
  analysis_data <- longitudinal_data[!is.na(get(biomarker)) & !is.na(liver_incident) & !is.na(follow_up_time)]
  if (nrow(analysis_data) < 50) next
  
  std_var <- paste0(biomarker,"_std")
  
  # Fit Cox model
  # Model includes same covariates as sensitivity analysis
  formula_str <- paste("Surv(follow_up_time, liver_incident) ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")
  cox_model <- coxph(as.formula(formula_str), data=analysis_data)
  
  # Check Proportional Hazards Assumption
  zph <- cox.zph(cox_model)
  
  # Extract global test result and result for the biomarker
  zph_table <- zph$table
  
  # Global test
  global_p <- zph_table["GLOBAL", "p"]
  
  # Biomarker test
  bm_p <- zph_table[std_var, "p"]
  
  cox_results <- rbind(cox_results, data.frame(
    biomarker = biomarker,
    biomarker_name = biomarker_names[biomarker],
    n_total = nrow(analysis_data),
    n_events = sum(analysis_data$liver_incident == 1),
    p_value_global = global_p,
    p_value_biomarker = bm_p,
    assumption_met = ifelse(bm_p > 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
}

# Save results
if (nrow(cox_results) > 0) {
  dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
  write.csv(cox_results, "results/charls/cox_ph_assumption_tests.csv", row.names=FALSE)
  print(cox_results)
}
