suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
data <- readRDS("data/processed/merged_individual_blood_2020.rds")
data <- data[inw1==1]
baseline_participants <- data[!is.na(r1livere)]
longitudinal_data <- baseline_participants[r1livere == 0]
longitudinal_data <- longitudinal_data[!is.na(nl_Blood_weight)]
longitudinal_data <- longitudinal_data[!is.na(r1agey)]
longitudinal_data <- longitudinal_data[!is.na(r1mbmi)]
longitudinal_data <- longitudinal_data[!is.na(r1smoken)]
longitudinal_data <- longitudinal_data[!is.na(r1drinkev)]
longitudinal_data[, liver_incident := ifelse((r2livere==1 & !is.na(r2livere)) | (r3livere==1 & !is.na(r3livere)) | (r4livere==1 & !is.na(r4livere)) | (r5livere==1 & !is.na(r5livere)), 1, 0)]
roc_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011","cti_mean","egdr_mean","tyg_mean","tyg_whtr_mean")
for (bm in roc_biomarkers) { std_col <- paste0(bm,"_std"); v <- longitudinal_data[[bm]]; idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident); mu <- mean(v[idx], na.rm=TRUE); sdv <- sd(v[idx], na.rm=TRUE); if (is.finite(sdv) && sdv != 0) longitudinal_data[, (std_col) := (get(bm)-mu)/sdv] else longitudinal_data[, (std_col) := NA_real_] }
biomarker_names_roc <- list("cti_2011"="CTI","egdr_2011"="eGDR","tyg_2011"="TyG","tyg_whtr_2011"="TyG-WHtR","cti_mean"="CTI (Mean)","egdr_mean"="eGDR (Mean)","tyg_mean"="TyG (Mean)","tyg_whtr_mean"="TyG-WHtR (Mean)")
roc_covariates <- c("r1agey","ragender","r1mbmi","raeducl","h1rural","marital","r1smoken","r1drinkev","r1chronic_group")
roc_cols <- c(roc_biomarkers, paste0(roc_biomarkers,"_std"), roc_covariates, "liver_incident")
roc_data_common <- longitudinal_data[complete.cases(longitudinal_data[, ..roc_cols])]

# Helper function for NRI and IDI
calculate_nri_idi <- function(outcome, prob_base, prob_new) {
  diff_prob <- prob_new - prob_base
  idi <- mean(diff_prob[outcome == 1]) - mean(diff_prob[outcome == 0])
  nri_cases <- mean(prob_new[outcome == 1] > prob_base[outcome == 1]) - mean(prob_new[outcome == 1] < prob_base[outcome == 1])
  nri_controls <- mean(prob_new[outcome == 0] < prob_base[outcome == 0]) - mean(prob_new[outcome == 0] > prob_base[outcome == 0])
  nri <- nri_cases + nri_controls
  return(list(nri = nri, idi = idi))
}

# Fit Base Model
base_model <- glm(liver_incident ~ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group, data=roc_data_common, family=binomial())
prob_base <- predict(base_model, type="response")

roc_results <- list()
auc_results <- data.frame(biomarker=character(), biomarker_name=character(), auc=numeric(), ci_lower=numeric(), ci_upper=numeric(), nri=numeric(), idi=numeric(), p_value=numeric(), n_total=numeric(), n_events=numeric(), stringsAsFactors=FALSE)
for (biomarker in roc_biomarkers) {
  current_data <- data.table::copy(roc_data_common)
  std_var <- paste0(biomarker,"_std")
  prediction_model <- glm(as.formula(paste("liver_incident ~", std_var, "+ r1agey + ragender + r1mbmi + raeducl + h1rural + marital + r1smoken + r1drinkev + r1chronic_group")), data=current_data, family=binomial())
  model_data <- prediction_model$model
  predicted_prob <- predict(prediction_model, type="response")
  actual_outcome <- model_data$liver_incident
  
  roc_obj <- roc(actual_outcome, predicted_prob, direction="auto", quiet=TRUE)
  auc_ci <- ci.auc(roc_obj, conf.level=0.95)
  roc_results[[biomarker]] <- roc_obj
  
  # NRI & IDI
  metrics <- calculate_nri_idi(actual_outcome, prob_base, predicted_prob)
  
  auc_results <- rbind(auc_results, data.frame(biomarker=biomarker, biomarker_name=biomarker_names_roc[[biomarker]], auc=as.numeric(auc_ci[2]), ci_lower=as.numeric(auc_ci[1]), ci_upper=as.numeric(auc_ci[3]), nri=metrics$nri, idi=metrics$idi, p_value=NA, n_total=length(actual_outcome), n_events=sum(actual_outcome==1, na.rm=TRUE), stringsAsFactors=FALSE))
}
if (length(roc_results) > 0) {
  roc_groups <- list("2011"=c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011"), "mean"=c("cti_mean","egdr_mean","tyg_mean","tyg_whtr_mean"))
  for (gname in names(roc_groups)) {
    gbiomarkers <- intersect(roc_groups[[gname]], names(roc_results))
    if (length(gbiomarkers)==0) next
    plot_data <- data.frame()
    for (biomarker in gbiomarkers) {
      roc_obj <- roc_results[[biomarker]]
      roc_coords <- coords(roc_obj, "all", ret=c("threshold","specificity","sensitivity"))
      temp_data <- data.frame(biomarker=biomarker, biomarker_name=biomarker_names_roc[[biomarker]], specificity=roc_coords$specificity, sensitivity=roc_coords$sensitivity, fpr=1-roc_coords$specificity)
      plot_data <- rbind(plot_data, temp_data)
    }
    p <- ggplot(plot_data, aes(x=fpr, y=sensitivity, color=biomarker_name)) + geom_line(size=1.2) + geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50") + scale_x_continuous(limits=c(0,1), expand=c(0,0)) + scale_y_continuous(limits=c(0,1), expand=c(0,0)) + labs(x="False positive rate (1 - Specificity)", y="Sensitivity", color="Prediction model") + theme_minimal() + theme(legend.position="bottom", legend.title=element_text(size=12), legend.text=element_text(size=10), axis.title=element_text(size=12), axis.text=element_text(size=10), panel.grid.minor=element_blank()) + guides(color=guide_legend(nrow=2, byrow=TRUE))
    if (nrow(auc_results) > 0) {
      auc_subset <- auc_results[auc_results$biomarker %in% gbiomarkers, ]
      if (nrow(auc_subset) > 0) {
        auc_labels <- paste0(auc_subset$biomarker_name, " (AUC: ", sprintf("%.3f", auc_subset$auc), ")")
        names(auc_labels) <- auc_subset$biomarker_name
        p <- p + scale_color_discrete(labels=auc_labels)
      }
    }
    dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
    outfile <- paste0("results/charls/roc_curves_", gname, ".tiff")
    ggsave(outfile, plot=p, width=10, height=8, dpi=300)
    print(p)
  }
}
if (nrow(auc_results) > 0) {
  cat("\n=== Prediction Model Performance Summary ===\n")
  for (i in 1:nrow(auc_results)) {
    cat(sprintf("%s:\n", auc_results$biomarker_name[i]))
    cat(sprintf("  Sample Size: %d (Events: %d)\n", auc_results$n_total[i], auc_results$n_events[i]))
    cat(sprintf("  AUC: %.3f (95%% CI: %.3f-%.3f)\n", auc_results$auc[i], auc_results$ci_lower[i], auc_results$ci_upper[i]))
    cat(sprintf("  NRI: %.3f\n", auc_results$nri[i]))
    cat(sprintf("  IDI: %.3f\n", auc_results$idi[i]))
    cat(sprintf("  P-value: %.3f\n\n", auc_results$p_value[i]))
  }
}

