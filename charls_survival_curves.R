suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
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
longitudinal_data[, follow_up_time := { time <- rep(0.5, .N); inc <- which(!is.na(liver_incident) & liver_incident == 1); for (i in inc) { if (!is.na(r2livere[i]) && r2livere[i] == 1) { time[i] <- 2 } else if (!is.na(r3livere[i]) && r3livere[i] == 1) { time[i] <- 4 } else if (!is.na(r4livere[i]) && r4livere[i] == 1) { time[i] <- 7 } else if (!is.na(r5livere[i]) && r5livere[i] == 1) { time[i] <- 9 } }; non <- which(is.na(liver_incident) | liver_incident == 0); for (i in non) { if (!is.na(r5livere[i])) { time[i] <- 9 } else if (!is.na(r4livere[i])) { time[i] <- 7 } else if (!is.na(r3livere[i])) { time[i] <- 4 } else if (!is.na(r2livere[i])) { time[i] <- 2 } }; time }]
longitudinal_data[follow_up_time <= 0, follow_up_time := 0.5]
main_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)")
kaplan_meier_plots <- list()
logrank_results <- data.frame(biomarker=character(), biomarker_name=character(), logrank_chisq=numeric(), logrank_p=numeric(), stringsAsFactors=FALSE)
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(longitudinal_data))) next
  analysis_data <- longitudinal_data[!is.na(get(biomarker)) & !is.na(liver_incident)]
  analysis_data[, time_to_event := follow_up_time]
  analysis_data[, event := liver_incident]
  analysis_data <- analysis_data[time_to_event > 0]
  if (nrow(analysis_data) >= 100) {
    qs <- quantile(analysis_data[[biomarker]], probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
    analysis_data[, biomarker_tertile := cut(get(biomarker), breaks=qs, labels=c("T1 (Low)", "T2 (Medium)", "T3 (High)"), include.lowest=TRUE)]
    surv_obj <- Surv(time=analysis_data$time_to_event, event=analysis_data$event)
    km_fit <- survfit(surv_obj ~ biomarker_tertile, data=analysis_data)
    lr <- survdiff(surv_obj ~ biomarker_tertile, data=analysis_data)
    lp <- 1 - pchisq(lr$chisq, df=length(lr$n) - 1)
    logrank_results <- rbind(logrank_results, data.frame(biomarker=biomarker, biomarker_name=biomarker_names[biomarker], logrank_chisq=lr$chisq, logrank_p=lp))
    km_plot <- ggsurvplot(km_fit, data=analysis_data, xlab="Time to Liver Disease (years)", ylab="Survival Probability", legend.title="Tertiles", legend.labs=c("T1 (Low)", "T2 (Medium)", "T3 (High)"), palette=c("#2E9FDF", "#00AFBB", "#E7B800"), risk.table=TRUE, risk.table.height=0.3, pval=TRUE, conf.int=FALSE, ggtheme=theme_minimal())
    dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
    ggsave(filename=paste0("results/charls/survival_curve_", biomarker, ".tiff"), plot=km_plot$plot, width=10, height=8, dpi=300)
    kaplan_meier_plots[[biomarker]] <- km_plot$plot + labs(title=biomarker_names[biomarker]) + theme(plot.title=element_text(size=14, face="bold"))
  }
}
write.csv(logrank_results, "results/charls/logrank_test_results.csv", row.names=FALSE)
if (length(kaplan_meier_plots) >= 4) {
  biomarker_order <- main_biomarkers[1:4]
  plots <- list()
  for (i in 1:4) {
    b <- biomarker_order[i]
    if (b %in% names(kaplan_meier_plots)) plots[[i]] <- kaplan_meier_plots[[b]]
  }
  combined_plot <- grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2, ncol=2)
  ggsave(filename="results/charls/survival_curves_combined_2x2.tiff", plot=combined_plot, width=16, height=12, dpi=300)
}

