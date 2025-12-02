suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))
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
wave_years <- c("r1"=2011,"r2"=2013,"r3"=2015,"r4"=2018,"r5"=2020)
longitudinal_data[, follow_up_time := { time <- rep(0.5, .N); incident_cases <- which(!is.na(liver_incident) & liver_incident == 1); for (i in incident_cases) { if (!is.na(r2livere[i]) && r2livere[i] == 1) { time[i] <- 2 } else if (!is.na(r3livere[i]) && r3livere[i] == 1) { time[i] <- 4 } else if (!is.na(r4livere[i]) && r4livere[i] == 1) { time[i] <- 7 } else if (!is.na(r5livere[i]) && r5livere[i] == 1) { time[i] <- 9 } }; non_incident <- which(is.na(liver_incident) | liver_incident == 0); for (i in non_incident) { if (!is.na(r5livere[i])) { time[i] <- 9 } else if (!is.na(r4livere[i])) { time[i] <- 7 } else if (!is.na(r3livere[i])) { time[i] <- 4 } else if (!is.na(r2livere[i])) { time[i] <- 2 } }; time }]
longitudinal_data[follow_up_time <= 0, follow_up_time := 0.5]
main_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)")
for (bm in main_biomarkers) { std <- paste0(bm,"_std"); v <- longitudinal_data[[bm]]; idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident); mu <- mean(v[idx], na.rm=TRUE); sdv <- sd(v[idx], na.rm=TRUE); if (is.finite(sdv) && sdv != 0) longitudinal_data[, (std) := (get(bm)-mu)/sdv] else longitudinal_data[, (std) := NA_real_] }
sensitivity_results <- data.frame()
sensitivity_data1 <- longitudinal_data[!is.na(liver_incident) & (is.na(r1hibpe) | r1hibpe == 0) & (is.na(r1diabe) | r1diabe == 0)]
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(sensitivity_data1))) next
  analysis_data <- sensitivity_data1[!is.na(get(biomarker))]
  if (nrow(analysis_data) < 50) next
  std_var <- paste0(biomarker,"_std")
  m <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev")), data=analysis_data, family=binomial())
  sm <- summary(m)
  or <- exp(sm$coefficients[std_var,"Estimate"])
  ci <- exp(confint(m)[std_var,])
  pv <- sm$coefficients[std_var,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Exclude other chronic diseases (Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(longitudinal_data))) next
  analysis_data <- longitudinal_data[!is.na(get(biomarker)) & !is.na(liver_incident)]
  if (nrow(analysis_data) < 50) next
  m <- glm(as.formula(paste("liver_incident ~", biomarker, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data, family=binomial())
  sm <- summary(m)
  or <- exp(sm$coefficients[biomarker,"Estimate"])
  ci <- exp(confint(m)[biomarker,])
  pv <- sm$coefficients[biomarker,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Continuous variable (Unstandardized, Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(longitudinal_data))) next
  analysis_data <- longitudinal_data[!is.na(get(biomarker)) & !is.na(liver_incident) & !is.na(follow_up_time)]
  if (nrow(analysis_data) < 50) next
  std_var <- paste0(biomarker,"_std")
  cox_model <- coxph(as.formula(paste("Surv(follow_up_time, liver_incident) ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data)
  cx <- summary(cox_model)
  hr <- cx$coefficients[std_var,"exp(coef)"]
  ci <- cx$conf.int[std_var,c("lower .95","upper .95")]
  pv <- cx$coefficients[std_var,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Cox regression (Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=hr, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
sensitivity_data_memory <- longitudinal_data[!is.na(liver_incident) & (is.na(r1memrye) | r1memrye == 0)]
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(sensitivity_data_memory))) next
  analysis_data <- sensitivity_data_memory[!is.na(get(biomarker))]
  if (nrow(analysis_data) < 50) next
  std_var <- paste0(biomarker,"_std")
  m <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data, family=binomial())
  sm <- summary(m)
  or <- exp(sm$coefficients[std_var,"Estimate"])
  ci <- exp(confint(m)[std_var,])
  pv <- sm$coefficients[std_var,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Exclude memory problems (Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
sensitivity_data_hbv <- longitudinal_data[!is.na(liver_incident) & (is.na(w1_hepatitis_b) | w1_hepatitis_b != 2)]
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(sensitivity_data_hbv))) next
  analysis_data <- sensitivity_data_hbv[!is.na(get(biomarker))]
  if (nrow(analysis_data) < 50) next
  std_var <- paste0(biomarker,"_std")
  m <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data, family=binomial())
  sm <- summary(m)
  or <- exp(sm$coefficients[std_var,"Estimate"])
  ci <- exp(confint(m)[std_var,])
  pv <- sm$coefficients[std_var,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Exclude Hepatitis B (Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
sensitivity_data_rxliver <- longitudinal_data[!is.na(liver_incident) & (is.na(r1rxliver_c) | r1rxliver_c == 0)]
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(sensitivity_data_rxliver))) next
  analysis_data <- sensitivity_data_rxliver[!is.na(get(biomarker))]
  if (nrow(analysis_data) < 50) next
  std_var <- paste0(biomarker,"_std")
  m <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data, family=binomial())
  sm <- summary(m)
  or <- exp(sm$coefficients[std_var,"Estimate"])
  ci <- exp(confint(m)[std_var,])
  pv <- sm$coefficients[std_var,"Pr(>|z|)"]
  sensitivity_results <- rbind(sensitivity_results, data.frame(analysis_type="Exclude liver-related medications (Longitudinal)", biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_cases=sum(analysis_data$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
}
if (nrow(sensitivity_results) > 0) {
  sensitivity_results$q_value <- ave(sensitivity_results$p_value, sensitivity_results$analysis_type, FUN=function(p) p.adjust(p, method="BH"))
  dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
  write.csv(sensitivity_results, "results/charls/sensitivity_analysis_longitudinal.csv", row.names=FALSE)
}

