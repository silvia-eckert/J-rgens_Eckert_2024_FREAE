# INFO ####
# This R script normalizes the raw UHPLC-QTOF-MS/MS data and adds metadata
# Author: Silvia Eckert
# Generated: 2024-05-13

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "dplyr", "tidyr"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Analysis_LCMS/R_tidy_LCMS_data.R") # Set root directory
setwd(here::here())
getwd()
library(dplyr)
library(tidyr)

# =========================================================================== #

# LOAD DATA ####
# metadata
meta <- read.table(here::here("Working_data",
                              "LCMS_data_metadata.txt"),
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = TRUE); head(meta); str(meta)
meta$Sample_ID <- as.character(meta$Sample_ID)

# bucket table
features <- read.table(here::here("Working_data",
                                "LCMS_data_raw.txt"),
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = TRUE); head(features); str(features)

# =========================================================================== #

# PREPARE DATA ####

# Step 1: Remove injections peaks (i.e. all buckets/features with RT <= 1.05 min)
features_tidy <- features %>% 
  dplyr::filter(RT > 1.05)

# Step 2: Divide by internal standard (row number 2821)
# Extract row with internal standard
# (Do not use the tidyverse here, as the row numbers are not retained!)
ISTD_row <- as.numeric(rownames(features_tidy[features_tidy$RT == 26.66 &
                                              features_tidy$m_z == 240.1033,])); ISTD_row # 2821
# Change unnecessary columns to character temporarily for convenience
features_tidy$RT <- as.character(features_tidy$RT)
features_tidy$m_z <- as.character(features_tidy$m_z)
# Divide by ISTD
features_tidy <- features_tidy %>% 
  mutate(across(where(is.numeric), ~ . / nth(.,n = ISTD_row)))
# ISTD row should have only 1s as values in columns starting with X...
features_tidy[features_tidy$RT == 26.66 & features_tidy$m_z == 240.1033,]

# Step 3: Subtract by blank sample column (and remove this column afterwards)
blank_column <- features_tidy %>%
  select(starts_with("X1_")) %>% 
  names() # X1_1207
# Subtract blank
features_tidy <- features_tidy %>% 
  mutate(across(where(is.numeric), ~ . - .data[[blank_column]])) %>% 
  dplyr::select(!X1_1207)

# Step 4: Replace negative values (from blank subtraction) with zeros
features_tidy[features_tidy < 0] <- 0

# Step 5: Rename samples
names(features_tidy) # before
features_tidy <- features_tidy %>% 
  dplyr::rename_with(~sub("X([0-9]+)_.*", "\\1", .))
features_tidy$RT <- as.numeric(features_tidy$RT)
features_tidy$m_z <- as.numeric(features_tidy$m_z)
features_tidy$Bucket_label <- as.character(features_tidy$Bucket_label)
names(features_tidy) # after

# Step 6: Transpose data
features_tidy <- features_tidy %>% 
  tidyr::pivot_longer(!c("Bucket_label", "RT", "m_z"),
                      names_to = "Sample_ID",
                      values_to = "peak_intensity")
  
# Step 7: Add metadata by sample name
features_clean <- left_join(meta, features_tidy, by = "Sample_ID") %>%
  dplyr::filter(Sample_type != "Blank")
# Check
head(features_clean); str(features_clean)

# Step 8: Add number for every level of RT
features_clean$bucket_id <- as.factor(features_clean$Bucket_label)
class(features_clean$bucket_id) # should be factor
levels(features_clean$bucket_id) <- c(1:length(levels(features_clean$bucket_id)))
features_clean$bucket_id # should be numeric factor

# Write results to file
write.table(features_clean,
            here::here("Working_data",
                       "LCMS_data_clean.txt"),
            row.names=F,
            quote=F,
            sep="\t")

# =========================T=H=E==E=N=D====================================== #
