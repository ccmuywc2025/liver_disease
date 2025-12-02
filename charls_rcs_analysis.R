suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rcssci))
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
main_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)")
for (bm in main_biomarkers) { std <- paste0(bm,"_std"); v <- longitudinal_data[[bm]]; idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident); mu <- mean(v[idx], na.rm=TRUE); sdv <- sd(v[idx], na.rm=TRUE); if (is.finite(sdv) && sdv != 0) longitudinal_data[, (std) := (get(bm)-mu)/sdv] else longitudinal_data[, (std) := NA_real_] }
df_analysis <- as.data.frame(longitudinal_data[!is.na(liver_incident)])
df_analysis$liver_outcome <- as.numeric(df_analysis$liver_incident)
biomarker_folders <- c("CTI","EGDR","TyG","TyG_WHtR")
names(biomarker_folders) <- main_biomarkers
for (folder_name in biomarker_folders) { dir.create(paste0("results/charls/rcs_plots/", folder_name), recursive=TRUE, showWarnings=FALSE); dir.create(paste0("results/charls/rcs_tables/", folder_name), recursive=TRUE, showWarnings=FALSE) }
for (biomarker in main_biomarkers) {
  folder_name <- biomarker_folders[biomarker]
  covariates <- c("r1agey","ragender","r1mbmi","raeducl","h1rural","marital","r1smoken","r1drinkev","r1chronic_group")
  available_covs <- covariates[covariates %in% names(df_analysis)]
  if (nrow(df_analysis) >= 100 && sum(df_analysis$liver_outcome == 1, na.rm=TRUE) >= 10) {
    res <- try(rcssci_logistic(data=df_analysis, y="liver_outcome", x=biomarker, covs=available_covs, knot=4, prob=0.1, filepath=paste0("results/charls/rcs_plots/", folder_name)), silent=TRUE)
    if (inherits(res, "try-error")) {
      for (knot_num in 3:5) {
        res <- try(rcssci_logistic(data=df_analysis, y="liver_outcome", x=biomarker, covs=available_covs, knot=knot_num, prob=0.1, filepath=paste0("results/charls/rcs_plots/", folder_name)), silent=TRUE)
        if (!inherits(res, "try-error")) { if (exists("result_table", res)) write.csv(res$result_table, paste0("results/charls/rcs_tables/", folder_name, "/rcs_logistic_", biomarker, "_knot", knot_num, ".csv"), row.names=FALSE); break }
      }
    } else { if (exists("result_table", res)) write.csv(res$result_table, paste0("results/charls/rcs_tables/", folder_name, "/rcs_logistic_", biomarker, ".csv"), row.names=FALSE) }
  }
}


