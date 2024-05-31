# INFO ####
# This R script summarizes the normalized UHPLC-QTOF-MS/MS data
# Author: Silvia Eckert
# Generated: 2024-05-13

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "dplyr", "tidyr", "vegan", "reshape2"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Analysis_LCMS/R_Tables_1_and_S3.R") # Set root directory
setwd(here::here())
getwd()
library(dplyr)
library(tidyr)
library(vegan)
library(reshape2)

# =========================================================================== #

# LOAD DATA ####
features <- read.table(here::here("Working_data",
                                  "LCMS_data_clean.txt"),
                       header = TRUE,
                       sep = "\t") %>%
  dplyr::select(c("Chemotype", "Treatment_ID", "Fraction_Methanol_ID", "Bucket_label",
                  "RT", "bucket_id", "m_z", "peak_intensity")) %>% 
  tidyr::pivot_wider(names_from = "Treatment_ID",
                     values_from = "peak_intensity",
                     values_fill = 0) %>% 
  dplyr::rename(Fraction = Fraction_Methanol_ID) %>% 
  dplyr::filter((Mock_treated + Low_dose + High_dose) != 0); features # 25,116 x 9

# =========================================================================== #

# RESTRUCTURE DATA ####
# Change to long-format structure
# NO FILTERING APPLIED!
features_long_unfiltered <- features %>%
  mutate(Fraction = factor(Fraction, levels = c("F10","F18","F30", "F100"))) %>%
  dplyr::select(c(Chemotype, Fraction, RT, bucket_id, Mock_treated, Low_dose, High_dose)) %>% 
  pivot_longer(cols = c(Mock_treated, Low_dose, High_dose),
               names_to = "Treatment",
               values_to = "Value") %>% 
  mutate(RT = factor(RT),
         bucket_id = factor(bucket_id),
         Treatment = factor(Treatment, levels = c("Mock_treated", "Low_dose", "High_dose")),
         Chemotype = as.factor(Chemotype)) %>% 
  dplyr::filter(Value > 0); features_long_unfiltered # 45,850 x 6

# Calculate difference compared to mock treatment for each 5AZA treatment level (low-dose and high-dose)
features_diff_CvsLD <- features %>% 
  dplyr::select(Chemotype, Fraction, bucket_id, Mock_treated, Low_dose) %>% 
  dplyr::filter((Mock_treated + Low_dose) > 0) %>% 
  mutate(diff = Low_dose - Mock_treated); features_diff_CvsLD # 20,293 x 6
features_diff_CvsHD <- features %>% 
  dplyr::select(Chemotype, Fraction, bucket_id, Mock_treated, High_dose) %>% 
  dplyr::filter((Mock_treated + High_dose) > 0) %>% 
  mutate(diff = High_dose - Mock_treated); features_diff_CvsHD # 21,034 x 6

# =========================================================================== #

# TABLE S3 ####
table_S3 <- features_long_unfiltered %>%
  group_by(Chemotype, Fraction, Treatment) %>%
  summarise(
    id_count_unique = n_distinct(bucket_id)) %>%
  ungroup() %>% 
  tidyr::pivot_wider(names_from = "Treatment",
                     values_from = "id_count_unique"); table_S3 # Table 1: 8 x 6

# Save table 1
if(!dir.exists(here("Tables"))){
  dir.create(here("Tables"))
}
write.table(table_S3,
            file = here("Tables",
                        "Table_S3.txt"),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

# =========================================================================== #

# TABLE 1 ####
# Step 1: Count higher/lower feature intensities per factor level combination
# C = mock treatment
# LD = low-dose treatment
# HD = high-dose treatment
# C vs LD
features_diff_counts_CvsLD <- features_diff_CvsLD %>%
  group_by(Chemotype, Fraction) %>%
  summarise(
    C_vs_LD.Higher = sum(diff > 0, na.rm = TRUE),
    C_vs_LD.Lower = sum(diff < 0, na.rm = TRUE),
    C_vs_LD.Sum = n()) %>% 
  tidyr::pivot_longer(cols = starts_with("C_"),
                      names_to = "Combination",
                      values_to = "Value") %>% 
  tidyr::separate(Combination,
                  into = c("Combination", "Features"),
                  sep = "\\.") %>% 
  tidyr::pivot_wider(names_from = Features,
                     values_from = Value) %>% 
  mutate(Higher_percent = Higher/Sum*100,
         Lower_percent = Lower/Sum*100); features_diff_counts_CvsLD

# C vs HD
features_diff_counts_CvsHD <- features_diff_CvsHD %>%
  group_by(Chemotype, Fraction) %>%
  summarise(
    C_vs_HD.Higher = sum(diff > 0, na.rm = TRUE),
    C_vs_HD.Lower = sum(diff < 0, na.rm = TRUE),
    C_vs_HD.Sum = n()) %>% 
  tidyr::pivot_longer(cols = starts_with("C_"),
                      names_to = "Combination",
                      values_to = "Value") %>% 
  tidyr::separate(Combination,
                  into = c("Combination", "Features"),
                  sep = "\\.") %>% 
  tidyr::pivot_wider(names_from = Features,
                     values_from = Value) %>% 
  mutate(Higher_percent = Higher/Sum*100,
         Lower_percent = Lower/Sum*100); features_diff_counts_CvsHD

# Combine count tables
features_diff_counts <- rbind(features_diff_counts_CvsLD,
                              features_diff_counts_CvsHD); features_diff_counts

# Step 2: Calculate Morisita-Horn index
# Prepare data
features_split <- features %>%
  dplyr::select(!c(Bucket_label, RT, bucket_id, m_z)) %>% 
  split(.$Chemotype)

features_Aketo <- features_split$Aketo %>%
  dplyr::select(!Chemotype) %>% 
  split(.$Fraction) %>% 
  purrr::map(~.x %>% select(-Fraction)) %>% 
  purrr::map(~.x %>% as.matrix()); features_Aketo

features_BThu <- features_split$BThu %>%
  dplyr::select(!Chemotype) %>% 
  split(.$Fraction) %>% 
  purrr::map(~.x %>% select(-Fraction)) %>% 
  purrr::map(~.x %>% as.matrix()); features_BThu

# Calculate Morisita-Horn similarity index
features_Aketo_mh <- data.frame()
for (fraction in names(features_Aketo)) {
  
  t_m <- t(features_Aketo[[fraction]])
  mh <- vegdist(t_m, method = "horn")
  
  mh_df <- melt(as.matrix(mh),
                varnames = c("row", "col")) %>% 
    dplyr::filter(!(row == col)) %>% 
    dplyr::filter(row_number() == 1:2) %>% 
    mutate(Chemotype = "Aketo",
           Fraction = fraction,
           .before = row) %>%
    mutate(row = dplyr::recode(row,
                               "Low_dose" = "LD",
                               "High_dose" = "HD"),
           col = dplyr::recode(col,
                               "Mock_treated" = "C")) %>% 
    mutate(Similarity = 1 - value) %>% 
    dplyr::select(!value) %>% 
    tidyr::unite(Combination,
                 c(col, row),
                 sep = "_vs_")
  features_Aketo_mh <- rbind(features_Aketo_mh, mh_df)
}; features_Aketo_mh

features_BThu_mh <- data.frame()
for (fraction in names(features_BThu)) {
  
  t_m <- t(features_BThu[[fraction]])
  mh <- vegdist(t_m, method = "horn")
  
  mh_df <- melt(as.matrix(mh),
                varnames = c("row", "col")) %>% 
    dplyr::filter(!(row == col)) %>% 
    dplyr::filter(row_number() == 1:2) %>% 
    mutate(Chemotype = "BThu",
           Fraction = fraction,
           .before = row) %>%
    mutate(row = dplyr::recode(row,
                               "Low_dose" = "LD",
                               "High_dose" = "HD"),
           col = dplyr::recode(col,
                               "Mock_treated" = "C")) %>%
    mutate(Similarity = 1 - value) %>% 
    dplyr::select(!value) %>% 
    tidyr::unite(Combination,
                 c(col, row),
                 sep = "_vs_")
  features_BThu_mh <- rbind(features_BThu_mh, mh_df)
}; features_BThu_mh

# Combine similarity tables
features_mh <- rbind(features_Aketo_mh,
                     features_BThu_mh); features_mh

# Step 3: Combine all tables and save
table_1 <- left_join(features_diff_counts %>% 
                       dplyr::select(Chemotype, Fraction, Combination, Higher_percent, Lower_percent, Sum) %>% 
                       dplyr::rename(Features = Sum),
                     features_mh) %>% 
  mutate(Combination = replace(Combination,
                               Combination == "C_vs_LD", "Mock_treated_vs_Low_dose"),
         Combination = replace(Combination,
                               Combination == "C_vs_HD", "Mock_treated_vs_High_dose")) %>% 
  dplyr::relocate(Chemotype, Combination, Fraction, Lower_percent, Higher_percent, Features, Similarity) %>%
  mutate(Fraction = factor(Fraction, levels = c("F10", "F18", "F30", "F100"))) %>% 
  mutate(Combination = sub('_vs_', ' vs ', Combination, perl = TRUE)) %>% 
  dplyr::arrange(Chemotype, rev(Combination), Fraction); table_1 # Table 1

# Save
if(!dir.exists(here("Tables"))){
  dir.create(here("Tables"))
}
write.table(table_1,
            file = here("Tables",
                        "Table_1.txt"),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

# =========================================================================== #

# COMPARE CHEMOTYPES ####
# Prepare data for overall comparison
chemotypes_pivoted <- features %>%
  tidyr::pivot_longer(cols = c("Mock_treated", "Low_dose", "High_dose"),
                      names_to = "Treatment",
                      values_to = "Intensity") %>% 
  tidyr::pivot_wider(names_from = Chemotype,
                     values_from = Intensity) %>% 
  mutate(Aketo = replace_na(Aketo, 0),
         BThu = replace_na(BThu, 0)) %>% 
  dplyr::select(c(Aketo, BThu)) %>% 
  as.matrix(); chemotypes_pivoted
# Prepare data for treatment-specific comparison of chemotypes
chemotypes_treatment_pivoted <- features %>%
  tidyr::pivot_longer(cols = c("Mock_treated", "Low_dose", "High_dose"),
                      names_to = "Treatment",
                      values_to = "Intensity") %>% 
  tidyr::pivot_wider(names_from = Chemotype,
                     values_from = Intensity) %>% 
  mutate(Aketo = replace_na(Aketo, 0),
         BThu = replace_na(BThu, 0)) %>% 
  dplyr::select(c(Treatment, Aketo, BThu)) %>% 
  split(.$Treatment) %>% 
  purrr::map(~.x %>% select(-Treatment)) %>% 
  purrr::map(~.x %>% as.matrix()); chemotypes_treatment_pivoted

# Calculate Morisita-Horn similarity index
# Overall
chemotypes_mh <- 1-vegdist(t(chemotypes_pivoted), method = "horn"); chemotypes_mh # 0.7903716
# Per treatment group
chemotypes_treatment_mh <- data.frame()
for (treatment in names(chemotypes_treatment_pivoted)) {
  subset <- chemotypes_treatment_pivoted[[treatment]]
  mh <- 1-vegdist(t(subset), method = "horn")
  result <- data.frame(Comparison = "Aketo_vs_BThu",
                       Treatment = treatment,
                       Similarity = mh)
  chemotypes_treatment_mh <- rbind(chemotypes_treatment_mh, result)
}
chemotypes_treatment_mh
#   Comparison    Treatment     Similarity
# 1 Aketo_vs_BThu High_dose     0.8401470
# 2 Aketo_vs_BThu Low_dose      0.7465612
# 3 Aketo_vs_BThu Mock_treated  0.8921892

# =========================T=H=E==E=N=D====================================== #
