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
insulin_related_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011","cti_mean","egdr_mean","tyg_mean","tyg_whtr_mean")
for (bm in insulin_related_biomarkers) { std_col <- paste0(bm,"_std"); v <- longitudinal_data[[bm]]; idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident); mu <- mean(v[idx], na.rm=TRUE); sdv <- sd(v[idx], na.rm=TRUE); if (is.finite(sdv) && sdv != 0) { longitudinal_data[, (std_col) := (get(bm)-mu)/sdv] } else { longitudinal_data[, (std_col) := NA_real_] } }
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)","cti_mean"="CTI (Mean)","egdr_mean"="eGDR (Mean)","tyg_mean"="TyG (Mean)","tyg_whtr_mean"="TyG-WHtR (Mean)")
longitudinal_results <- data.frame(biomarker=character(), biomarker_name=character(), n_total=numeric(), n_events=numeric(), or_crude=numeric(), ci_lower_crude=numeric(), ci_upper_crude=numeric(), p_crude=numeric(), or_adjusted=numeric(), ci_lower_adjusted=numeric(), ci_upper_adjusted=numeric(), p_adjusted=numeric(), stringsAsFactors=FALSE)
for (biomarker in insulin_related_biomarkers) {
  analysis_data <- longitudinal_data[!is.na(get(biomarker)) & !is.na(liver_incident)]
  if (nrow(analysis_data) < 50 || sum(analysis_data$liver_incident==1, na.rm=TRUE) < 10) next
  std_var <- paste0(biomarker,"_std")
  model_crude <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender")), data=analysis_data, family=binomial())
  model_adjusted <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=analysis_data, family=binomial())
  crude_summary <- summary(model_crude)
  adjusted_summary <- summary(model_adjusted)
  or_crude <- exp(crude_summary$coefficients[std_var,"Estimate"])
  ci_crude <- exp(confint(model_crude)[std_var,])
  p_crude <- crude_summary$coefficients[std_var,"Pr(>|z|)"]
  or_adjusted <- exp(adjusted_summary$coefficients[std_var,"Estimate"])
  ci_adjusted <- exp(confint(model_adjusted)[std_var,])
  p_adjusted <- adjusted_summary$coefficients[std_var,"Pr(>|z|)"]
  longitudinal_results <- rbind(longitudinal_results, data.frame(biomarker=biomarker, biomarker_name=biomarker_names[biomarker], n_total=nrow(analysis_data), n_events=sum(analysis_data$liver_incident==1, na.rm=TRUE), or_crude=or_crude, ci_lower_crude=ci_crude[1], ci_upper_crude=ci_crude[2], p_crude=p_crude, or_adjusted=or_adjusted, ci_lower_adjusted=ci_adjusted[1], ci_upper_adjusted=ci_adjusted[2], p_adjusted=p_adjusted, stringsAsFactors=FALSE))
}
if (nrow(longitudinal_results) > 0) {
  primary_bms <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
  secondary_bms <- c("cti_mean","egdr_mean","tyg_mean","tyg_whtr_mean")
  longitudinal_results$q_value <- NA_real_
  idx1 <- longitudinal_results$biomarker %in% primary_bms
  idx2 <- longitudinal_results$biomarker %in% secondary_bms
  if (any(idx1, na.rm=TRUE)) longitudinal_results$q_value[idx1] <- p.adjust(longitudinal_results$p_adjusted[idx1], method="BH")
  if (any(idx2, na.rm=TRUE)) longitudinal_results$q_value[idx2] <- p.adjust(longitudinal_results$p_adjusted[idx2], method="BH")
  dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
  write.csv(longitudinal_results, "results/charls/novel_biomarkers_longitudinal_results.csv", row.names=FALSE)
}


