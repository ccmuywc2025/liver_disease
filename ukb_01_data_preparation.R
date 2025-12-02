# UKB External Validation Data Preparation Script
# Function: Read raw data, perform cleaning, exclusion, variable derivation, and standardization
# Input: Raw CSV data files
# Output: Processed RDS data file (data/processed/ukb_validation_dataset.rds)


source("R/code_release/ukb/ukb_00_functions.R")
suppressPackageStartupMessages(library(dplyr))

# 1) Input and Output Paths
participant_path <- "data/ukb数据/肝病_participant.csv"
output_rds_path <- "data/processed/ukb_validation_dataset.rds"
output_csv_path <- "data/processed/ukb_validation_summary.csv"

# 2) Read Participant Data (UTF-8+BOM compatible)
if (!file.exists(participant_path)) {
  stop(paste("File does not exist:", participant_path))
}
df <- read.csv(participant_path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)

supplementary_data_path <- "data/ukb数据/肝病结局来源和mri补充_participant.csv"
if (file.exists(supplementary_data_path)) {
  supp <- read.csv(supplementary_data_path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, supp, by = "eid")
}

# 3) Convert UKB Date Columns and Establish Project-Specific Naming
date_cols <- intersect(
  c(
    "p53_i0",   # Baseline visit date
    "p131286",  # I10 First reported date
    "p131294",  # I15 First reported date
    "p131658",  # K70 First reported date
    "p131660",  # K71 First reported date
    "p131662",  # K72 First reported date
    "p131664",  # K73 First reported date
    "p131666",  # K74 First reported date
    "p131670",  # K76 First reported date
    "p131668"   # K75 First reported date
  ),
  names(df)
)
df[date_cols] <- lapply(df[date_cols], as.Date)

dlc_vec <- do.call(c, df[date_cols])
data_lock_date <- if (all(is.na(dlc_vec))) as.Date(NA) else suppressWarnings(max(dlc_vec, na.rm = TRUE))

# Baseline Liver Disease Exclusion: If any liver disease first reported date <= baseline date, consider as baseline liver disease
event_cols <- intersect(c("p131658","p131660","p131662","p131664","p131666","p131670","p131668"), names(df))
df <- df %>%
  mutate(
    baseline_cld = {
      if (length(event_cols) == 0) {
        as.integer(rep(NA, n()))
      } else {
        flags <- as.data.frame(across(all_of(event_cols), ~ !is.na(.x) & !is.na(p53_i0) & (.x <= p53_i0)))
        cnt <- rowSums(flags)
        ifelse(is.na(p53_i0), NA_integer_, ifelse(cnt > 0, 1L, 0L))
      }
    }
  )
excluded_baseline_cld_n <- sum(df$baseline_cld %in% 1L, na.rm = TRUE)
cat(sprintf("Number of participants excluded with baseline liver disease: %d\n", excluded_baseline_cld_n))
df <- df %>% dplyr::filter(!(baseline_cld %in% 1L))

# 3.5) Exclude Liver Cancer Patients from Supplementary Data
supplementary_path <- "data/ukb数据/肝病队列癌症补充_participant.csv"
if (file.exists(supplementary_path)) {
  cat("Processing liver cancer patients in supplementary data...\n")
  supp_data <- read.csv(supplementary_path, stringsAsFactors = FALSE)
  cancer_fields <- grep("^(p20001|p20007)", names(supp_data), value = TRUE)
  
  if (length(cancer_fields) > 0) {
    liver_cancer_patterns <- c("liver/hepatocellular cancer", "hepatocellular", "liver cancer", "C22")
    has_liver_cancer <- rep(FALSE, nrow(supp_data))
    for (field in cancer_fields) {
      if (field %in% names(supp_data)) {
        for (pattern in liver_cancer_patterns) {
          has_liver_cancer <- has_liver_cancer | 
            grepl(pattern, supp_data[[field]], ignore.case = TRUE, fixed = FALSE)
        }
      }
    }
    excluded_eids <- supp_data$eid[has_liver_cancer]
    n_excluded_liver_cancer <- length(excluded_eids)
    cat(sprintf("Found liver cancer patients in supplementary data: %d\n", n_excluded_liver_cancer))
    if (n_excluded_liver_cancer > 0 && "eid" %in% names(df)) {
      df <- df %>% dplyr::filter(!eid %in% excluded_eids)
      cat(sprintf("Excluded liver cancer patients from supplementary data, remaining participants: %d\n", nrow(df)))
    }
  }
}

# 4) Unit Conversion (TG/GLU -> mg/dL)
df <- df %>%
  mutate(
    tg_mgdl = as.numeric(p30870_i0) * 88.57,
    glu_mgdl = as.numeric(p30740_i0) * 18,
    hba1c_pct = (as.numeric(p30750_i0) / 10.929) + 2.15
  )

# 5) Baseline Hypertension
df <- df %>%
  mutate(
    baseline_htn = case_when(
      is.na(p53_i0) ~ as.integer(NA),
      (!is.na(p131286) & p131286 <= p53_i0) | (!is.na(p131294) & p131294 <= p53_i0) ~ 1L,
      TRUE ~ 0L
    )
  )

# 6) Calculate Four Indices (CTI, eGDR, TyG, WHtR, TyG-WHtR)
df <- df %>%
  mutate(
    cti_ukb = 0.412 * log(p30710_i0) + (log(tg_mgdl * glu_mgdl))/2,
    egdr_ukb = 21.158 - (0.09 * p48_i0) - (3.407 * baseline_htn) - (0.551 * hba1c_pct),
    whtr_ukb = p48_i0 / p50_i0,
    tyg_ukb = log((tg_mgdl * glu_mgdl) / 2),
    tyg_whtr_ukb = tyg_ukb * whtr_ukb
  )

# 7) Outcome Processing
after_cols <- paste0(event_cols, "_after")
df <- df %>%
  mutate(
    across(
      all_of(event_cols),
      ~ if_else(!is.na(.x) & !is.na(p53_i0) & .x >= p53_i0, .x, as.Date(NA)),
      .names = "{.col}_after"
    )
  ) %>%
  mutate(
    cld_event_date = {
      if (length(after_cols) == 0) {
        as.Date(rep(NA, n()))
      } else {
        cols <- across(all_of(after_cols))
        out <- suppressWarnings(do.call(pmin, c(cols, na.rm = TRUE)))
        if (!inherits(out, "Date")) {
          out <- as.Date(out, origin = "1970-01-01")
        }
        na_all <- rowSums(!is.na(cols)) == 0L
        out[na_all] <- as.Date(NA)
        out
      }
    }
  ) %>%
  mutate(
    incident_cld = if_else(is.na(cld_event_date), 0L, 1L)
  )

df <- df %>%
  mutate(
    follow_up_end_date = if_else(incident_cld == 1L & !is.na(cld_event_date), cld_event_date, data_lock_date),
    follow_up_time_days = if_else(is.na(p53_i0) | is.na(follow_up_end_date), as.numeric(NA), pmax(0, as.numeric(difftime(follow_up_end_date, p53_i0, units = "days")))),
    follow_up_time_years = if_else(is.na(follow_up_time_days), as.numeric(NA), follow_up_time_days / 365.25)
  )

df <- df %>%
  mutate(
    across(
      all_of(after_cols),
      ~ if_else(is.na(.x), 0L, 1L),
      .names = "{.col}_event"
    )
  ) %>%
  mutate(
    across(
      all_of(after_cols),
      ~ if_else(!is.na(.x), .x, data_lock_date),
      .names = "{.col}_end"
    )
  ) %>%
  mutate(
    across(
      all_of(paste0(after_cols, "_end")),
      ~ if_else(is.na(p53_i0) | is.na(.x), as.numeric(NA), pmax(0, as.numeric(difftime(.x, p53_i0, units = "days")))),
      .names = "{.col}_time_days"
    )
  ) %>%
  mutate(
    across(
      all_of(paste0(after_cols, "_end")),
      ~ as.numeric(if_else(is.na(p53_i0) | is.na(.x), NA_real_, pmax(0, as.numeric(difftime(.x, p53_i0, units = "days"))) / 365.25)),
      .names = "{.col}_time_years"
    )
  )

# Covariate Processing
df <- df %>%
  mutate(
    p20116_i0_3cat = case_when(
      p20116_i0 == "Current" ~ "Current",
      p20116_i0 == "Previous" ~ "Previous",
      p20116_i0 == "Never"  ~ "Never", 
      TRUE ~ NA_character_
    ),
    p1558_i0_3cat = case_when(
      p1558_i0 %in% c("Daily or almost daily", "Three or four times a week") ~ "High",
      p1558_i0 %in% c("Once or twice a week", "One to three times a month") ~ "Medium",
      p1558_i0 %in% c("Special occasions only", "Never") ~ "Low",
      p1558_i0 %in% c("Prefer not to answer", "") | is.na(p1558_i0) ~ NA_character_, 
      TRUE ~ NA_character_
    ),
    p6138_i0_3cat = case_when(
      p6138_i0 %in% c("None of the above", "CSEs or equivalent") |
        grepl("CSEs or equivalent", p6138_i0) & !grepl("College or University degree", p6138_i0) ~ "Lower_Education",
      grepl("College or University degree", p6138_i0) ~ "Higher_Education",
      p6138_i0 %in% c("O levels/GCSEs or equivalent", "A levels/AS levels or equivalent", "NVQ or HND or HNC or equivalent", "Other professional qualifications eg: nursing, teaching") |
        grepl("O levels/GCSEs or equivalent|A levels/AS levels or equivalent|NVQ or HND or HNC or equivalent", p6138_i0) & 
        !grepl("College or University degree", p6138_i0) ~ "Intermediate_Education",
      p6138_i0 %in% c("Prefer not to answer", "") | is.na(p6138_i0) ~ NA_character_,
      TRUE ~ NA_character_
    )
  )
df$p6138_i0_3cat <- factor(df$p6138_i0_3cat, levels = c("Lower_Education", "Intermediate_Education", "Higher_Education"))

df <- df %>%
  mutate(
    p20118_i0_urban_rural = case_when(
      p20118_i0 %in% c("England/Wales - Urban - less sparse", "England/Wales - Urban - sparse", "Scotland - Large Urban Area", "Scotland - Other Urban Area") ~ "Urban",
      p20118_i0 %in% c("England/Wales - Town and Fringe - less sparse", "England/Wales - Town and Fringe - sparse", "England/Wales - Village - less sparse", "England/Wales - Village - sparse", "England/Wales - Hamlet and Isolated Dwelling - less sparse", "England/Wales - Hamlet and Isolated dwelling - sparse", "Scotland - Accessible Small Town", "Scotland - Accessible Rural", "Scotland - Remote Small Town", "Scotland - Remote Rural", "Scotland - Very Remote Rural") ~ "Rural",
      p20118_i0 %in% c("", "Postcode not linkable") | is.na(p20118_i0) ~ NA_character_,
      TRUE ~ NA_character_
    )
  )
df$p20118_i0_urban_rural <- factor(df$p20118_i0_urban_rural)

# Winsorize Processing
continuous_vars <- c("cti_ukb", "egdr_ukb", "tyg_ukb", "tyg_whtr_ukb")
cat("\n1. Extreme Value Processing:\n")
for (var in continuous_vars) {
  if (var %in% names(df)) {
    before_mean <- mean(df[[var]], na.rm = TRUE)
    df[[var]] <- winsorize_variable(df[[var]])
    after_mean <- mean(df[[var]], na.rm = TRUE)
    cat(sprintf("  %s: Mean changed from %.3f to %.3f\n", var, before_mean, after_mean))
  }
}

# Standardization
cat("\n2. Standardized Variable Creation:\n")
standardized_vars <- c()
for (var in continuous_vars) {
  if (var %in% names(df)) {
    std_var_name <- paste0(var, "_std")
    df[[std_var_name]] <- standardize_variable(df[[var]])
    standardized_vars <- c(standardized_vars, std_var_name)
    cat(sprintf("  %s standardization completed\n", var))
  }
}

# Missing Value Handling
available_covariates <- c(
  "p21022", "p31", "p21001_i0", "p20116_i0_3cat", "p1558_i0_3cat", "p6138_i0_3cat", "p20118_i0_urban_rural"
)
available_covariates <- available_covariates[available_covariates %in% names(df)]

cat("\n3. Missing Value Exclusion:\n")
df <- df %>%
  dplyr::filter(complete.cases(across(all_of(c(continuous_vars, available_covariates)))))

# 8) Save RDS
dir.create(dirname(output_rds_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(df, file = output_rds_path)
cat(sprintf("Data preparation completed, saved to %s\n", output_rds_path))
