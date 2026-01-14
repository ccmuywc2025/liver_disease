suppressPackageStartupMessages(library(data.table))
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
longitudinal_data[, liver_incident := ifelse((r2livere == 1 & !is.na(r2livere)) | (r3livere == 1 & !is.na(r3livere)) | (r4livere == 1 & !is.na(r4livere)) | (r5livere == 1 & !is.na(r5livere)), 1, 0)]
main_biomarkers <- c("cti_2011","egdr_2011","tyg_2011","tyg_whtr_2011")
for (bm in main_biomarkers) { std <- paste0(bm,"_std"); v <- longitudinal_data[[bm]]; idx <- !is.na(v) & !is.na(longitudinal_data$liver_incident); mu <- mean(v[idx], na.rm=TRUE); sdv <- sd(v[idx], na.rm=TRUE); if (is.finite(sdv) && sdv != 0) longitudinal_data[, (std) := (get(bm)-mu)/sdv] else longitudinal_data[, (std) := NA_real_] }
biomarker_names <- c("cti_2011"="CTI(2011)","egdr_2011"="eGDR(2011)","tyg_2011"="TyG(2011)","tyg_whtr_2011"="TyG-WHtR(2011)")
subgroup_vars <- list(
  "Gender"=list(var="ragender",levels=c("Male"=1,"Female"=2)),
  "Age Group"=list(var="age_group",levels=c("<60 years"=0,">=60 years"=1)),
  "BMI Group"=list(var="bmi_group",levels=c("Normal/Lean"=0,"Overweight/Obese"=1)),
  "Residence"=list(var="h1rural",levels=c("Urban"=0,"Rural"=1))
)
subgroup_data <- longitudinal_data[!is.na(liver_incident)]
subgroup_data[, age_group := ifelse(r1agey >= 60, 1, 0)]
subgroup_data[, bmi_group := ifelse(r1mbmi >= 24, 1, 0)]
subgroup_results <- data.frame()
for (biomarker in main_biomarkers) {
  if (!(biomarker %in% names(subgroup_data))) next
  for (subgroup_name in names(subgroup_vars)) {
    subgroup_var <- subgroup_vars[[subgroup_name]]$var
    if (!(subgroup_var %in% names(subgroup_data))) next
    for (level_name in names(subgroup_vars[[subgroup_name]]$levels)) {
      level_value <- subgroup_vars[[subgroup_name]]$levels[[level_name]]
      subset_dt <- subgroup_data[get(subgroup_var) == level_value & !is.na(get(biomarker))]
      if (nrow(subset_dt) < 50) next
      std_var <- paste0(biomarker,"_std")
      base_covariates <- c("r1agey","ragender","r1mbmi","raeducl","h1rural","marital","r1smoken","r1drinkev","r1chronic_group")
      if (subgroup_var == "ragender") covariates <- base_covariates[!base_covariates %in% c("ragender")] else if (subgroup_var == "age_group") covariates <- base_covariates[!base_covariates %in% c("r1agey")] else if (subgroup_var == "bmi_group") covariates <- base_covariates[!base_covariates %in% c("r1mbmi")] else if (subgroup_var == "h1rural") covariates <- base_covariates[!base_covariates %in% c("h1rural")] else covariates <- base_covariates
      fml <- as.formula(paste("liver_incident ~", std_var, "+", paste(covariates, collapse=" + ")))
      md <- glm(fml, data=subset_dt, family=binomial())
      sm <- summary(md)
      or <- exp(sm$coefficients[std_var,"Estimate"])
      ci <- exp(confint(md)[std_var,])
      pv <- sm$coefficients[std_var,"Pr(>|z|)"]
      subgroup_results <- rbind(subgroup_results, data.frame(biomarker=biomarker, biomarker_name=biomarker_names[biomarker], subgroup=subgroup_name, subgroup_level=level_name, n_total=nrow(subset_dt), n_cases=sum(subset_dt$liver_incident==1), or=or, ci_lower=ci[1], ci_upper=ci[2], p_value=pv, stringsAsFactors=FALSE))
    }
  }
}
if (nrow(subgroup_results) > 0) {
  subgroup_results$q_value <- ave(subgroup_results$p_value, subgroup_results$biomarker, FUN=function(p) p.adjust(p, method="BH"))
  dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
  write.csv(subgroup_results, "results/charls/subgroup_analysis_longitudinal.csv", row.names=FALSE)
}
interaction_df <- NULL
if (file.exists("results/charls/interaction_test_results.csv")) interaction_df <- tryCatch(read.csv("results/charls/interaction_test_results.csv", stringsAsFactors=FALSE), error=function(e) NULL)
biomarker_names_english <- c("cti_2011"="CTI","egdr_2011"="eGDR","tyg_2011"="TyG","tyg_whtr_2011"="TyG-WHtR")
subgroup_category_map <- data.frame(
  subgroup=c("Age Group","Gender","BMI Group","Residence"),
  subgroup_english=c("Age Group","Gender","BMI Group","Residence"),
  category=c("Age Group","Gender","BMI Group","Residence"),
  order=c(1,2,3,5)
)
subgroup_level_map <- data.frame(
  subgroup_level=c("<60 years",">=60 years","Male","Female","Normal/Lean","Overweight/Obese","Urban","Rural"),
  level_english=c("<60 years",">=60 years","Male","Female","Normal/Lean","Overweight/Obese","Urban","Rural")
)
plot_data <- subgroup_results %>% left_join(subgroup_category_map, by="subgroup") %>% left_join(subgroup_level_map, by="subgroup_level") %>% mutate(subgroup_english=ifelse(is.na(subgroup_english), subgroup, subgroup_english), level_english=ifelse(is.na(level_english), subgroup_level, level_english), category=ifelse(is.na(category), subgroup_english, category), order=ifelse(is.na(order), 99, order), subgroup_label=paste(subgroup_english, level_english, sep=": "), biomarker_label=biomarker_names_english[biomarker], significant=ifelse(p_value<0.05, "Significant", "Non-significant"), ci_text=sprintf("%.2f (%.2f-%.2f)", or, ci_lower, ci_upper), category_order=paste0(sprintf("%02d", order), "_", category), subgroup_order=paste0(category_order, "_", subgroup_english, "_", level_english)) %>% arrange(category_order, subgroup_order)
interaction_labels <- NULL
if (!is.null(interaction_df) && nrow(interaction_df) > 0) {
  interaction_labels <- plot_data %>% left_join(interaction_df %>% select(biomarker_name, subgroup, interaction_p_lrt) %>% rename(interaction_p=interaction_p_lrt), by=c("biomarker_name","subgroup")) %>% filter(!is.na(interaction_p)) %>% group_by(biomarker_label, subgroup_english) %>% summarise(interaction_p=first(interaction_p), max_or=max(or, na.rm=TRUE), max_ci_upper=max(ci_upper, na.rm=TRUE), subgroup_label=first(subgroup_label), .groups="drop") %>% mutate(interaction_text=ifelse(interaction_p<0.05, sprintf("P-int=%.3f*", interaction_p), sprintf("P-int=%.3f", interaction_p)), label_x=max_ci_upper*1.2, label_y=subgroup_label)
}
p1 <- ggplot(plot_data, aes(x=or, y=reorder(subgroup_label, -as.numeric(factor(subgroup_order))))) + geom_vline(xintercept=1, linetype="dashed", color="red", alpha=0.7) + geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper, color=significant), height=0.15, size=0.8) + geom_point(aes(color=significant, size=n_total), alpha=0.8) + geom_text(aes(label=ci_text), hjust=-0.1, size=2.5) + scale_color_manual(values=c("Significant"="#E31A1C","Non-significant"="#1F78B4")) + scale_size_continuous(range=c(2,4), name="Sample Size") + scale_x_log10(breaks=c(0.5,0.7,1.0,1.5,2.0,3.0), labels=c("0.5","0.7","1.0","1.5","2.0","3.0")) + facet_grid(category ~ biomarker_name, scales="free_y", space="free_y") + theme_minimal() + theme(legend.position="bottom", panel.grid.minor=element_blank(), strip.text.x=element_text(face="bold", size=10), strip.text.y=element_text(face="bold", size=9, angle=0))
if (!is.null(interaction_labels) && nrow(interaction_labels) > 0) p1 <- p1 + geom_text(data=interaction_labels, aes(x=label_x, y=subgroup_label, label=interaction_text), hjust=0, size=2.2, color="darkblue", fontface="italic", inherit.aes=FALSE)
ggsave("results/charls/subgroup_analysis_facet_plot_improved.tiff", plot=p1, width=16, height=12, dpi=300)
p2 <- ggplot(plot_data, aes(x=or, y=reorder(subgroup_label, -as.numeric(factor(subgroup_order))))) + geom_vline(xintercept=1, linetype="dashed", color="red", alpha=0.7) + geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper, color=significant), height=0.15, size=0.8) + geom_point(aes(color=significant), alpha=0.8, size=3) + scale_color_manual(values=c("Significant"="#E31A1C","Non-significant"="#1F78B4")) + scale_x_log10(breaks=c(0.5,0.7,1.0,1.5,2.0,3.0), labels=c("0.5","0.7","1.0","1.5","2.0","3.0")) + facet_wrap(~ biomarker_label, scales="free_y", ncol=2, nrow=2) + theme_minimal() + theme(legend.position="bottom", strip.text=element_text(face="bold", size=10), panel.grid.minor=element_blank(), panel.grid.major.y=element_line(color="grey90", size=0.3), panel.grid.major.x=element_line(color="grey90", size=0.3))
if (!is.null(interaction_labels) && nrow(interaction_labels) > 0) p2 <- p2 + geom_text(data=interaction_labels, aes(x=label_x, y=label_y, label=interaction_text), hjust=1, vjust=1, size=3, color="darkblue", fontface="italic", inherit.aes=FALSE)
ggsave("results/charls/subgroup_analysis_facet_plot_2x2.tiff", plot=p2, width=12, height=10, dpi=300)
summary_table <- plot_data %>% select(biomarker_name, subgroup, subgroup_level, n_total, n_cases, or, ci_lower, ci_upper, p_value) %>% mutate(or_ci=sprintf("%.2f (%.2f-%.2f)", or, ci_lower, ci_upper), p_formatted=ifelse(p_value<0.001, "<0.001", sprintf("%.3f", p_value))) %>% select(biomarker_name, subgroup, subgroup_level, n_total, n_cases, or_ci, p_formatted) %>% rename("Biomarker"=biomarker_name, "Subgroup"=subgroup, "Subgroup Level"=subgroup_level, "Total Sample"=n_total, "Cases"=n_cases, "OR (95%CI)"=or_ci, "P-value"=p_formatted)
dir.create("results/charls", showWarnings=FALSE, recursive=TRUE)
write.csv(summary_table, "results/charls/subgroup_analysis_summary_table.csv", row.names=FALSE)

