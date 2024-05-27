# INFO ####
# This R script plots the normalized UHPLC-QTOF-MS/MS data as a heatmap
# Author: Silvia Eckert
# Generated: 2024-05-13

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "ggplot2", "dplyr", "tidyr", "ggpubr", "cowplot"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Analysis_LCMS/R_Figure_2.R") # Set root directory
setwd(here::here())
getwd()
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)

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
features_long <- features %>% 
  filter(Mock_treated > 0.5, Low_dose > 0.5, High_dose > 0.5) %>% # 185 × 9 unique features
  mutate(Fraction = factor(Fraction, levels = c("F10","F18","F30", "F100"))) %>%
  dplyr::select(c(Chemotype, Fraction, RT, bucket_id, Mock_treated, Low_dose, High_dose)) %>% 
  pivot_longer(cols = c(Mock_treated, Low_dose, High_dose),
               names_to = "Treatment",
               values_to = "Intensity") %>% 
  mutate(Intensity_scaled = scale(Intensity)[,1],
         bucket_id = factor(bucket_id),
         Treatment = factor(Treatment, levels = c("Mock_treated", "Low_dose", "High_dose")),
         Chemotype = as.factor(Chemotype)); features_long # 555 (185 features * 3 treatment groups) x 7

# =========================================================================== #

# SUBSET DATA ####
# find intersection of both chemotypes per fraction
# Fraction F10
Aketo_F10_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "Aketo" &
                                                          features_long$Fraction == "F10"])
BThu_F10_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "BThu" &
                                                         features_long$Fraction == "F10"])
F10_id_intersect <- as.numeric(intersect(Aketo_F10_id_unique,
                                         BThu_F10_id_unique)); F10_id_intersect; length(F10_id_intersect) # 12
# Fraction F18
Aketo_F18_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "Aketo" &
                                                          features_long$Fraction == "F18"])
BThu_F18_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "BThu" &
                                                         features_long$Fraction == "F18"])
F18_id_intersect <- as.numeric(intersect(Aketo_F18_id_unique,
                                         BThu_F18_id_unique)); F18_id_intersect; length(F18_id_intersect) # 24
# Fraction F30
Aketo_F30_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "Aketo" &
                                                          features_long$Fraction == "F30"])
BThu_F30_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "BThu" &
                                                         features_long$Fraction == "F30"])
F30_id_intersect <- as.numeric(intersect(Aketo_F30_id_unique,
                                         BThu_F30_id_unique)); F30_id_intersect; length(F30_id_intersect) # 8
# Fraction F100
Aketo_F100_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "Aketo" &
                                                           features_long$Fraction == "F100"])
BThu_F100_id_unique <- as.character(features_long$bucket_id[features_long$Chemotype == "BThu" &
                                                          features_long$Fraction == "F100"])
F100_id_intersect <- as.numeric(intersect(Aketo_F100_id_unique,
                                          BThu_F100_id_unique)); F100_id_intersect; length(F100_id_intersect) # 15

# Subset to intersection for each fraction
features_F10_subset <- features_long[features_long$bucket_id %in% F10_id_intersect &
                               features_long$Fraction == "F10",]; features_F10_subset
features_F18_subset <- features_long[features_long$bucket_id %in% F18_id_intersect &
                               features_long$Fraction == "F18",]; features_F18_subset
features_F30_subset <- features_long[features_long$bucket_id %in% F30_id_intersect &
                               features_long$Fraction == "F30",]; features_F30_subset
features_F100_subset <- features_long[features_long$bucket_id %in% F100_id_intersect &
                               features_long$Fraction == "F100",]; features_F100_subset
# Combine
features_subset <- rbind(features_F10_subset, features_F18_subset,
                     features_F30_subset, features_F100_subset); features_subset # 354 x 7
# Check if identical
identical(features_subset$bucket_id[features_subset$Chemotype == "Aketo"],
          features_subset$bucket_id[features_subset$Chemotype == "BThu"]) # should be TRUE
  
# =========================================================================== #

# PLOT DATA ####
# Total number of features in heatmap
length(unique(features_subset$bucket_id)) # 54

# Heatmap for chemotype Aketo
p1 <- ggplot(features_subset %>%
               dplyr::filter(Chemotype == "Aketo") %>% 
               mutate(Treatment = forcats::fct_recode(Treatment,
                                                      "Mock-treated" = "Mock_treated",
                                                      "Low-dose" = "Low_dose",
                                                      "High-dose" = "High_dose")),
             aes(y = bucket_id, x = Treatment, fill = Intensity_scaled)) + 
  geom_tile() +
  labs(y = "Feature ID",
       fill = "Standardized and scaled feature intensity") +
  facet_grid(Fraction ~ Chemotype, scales = "free") +
  scale_fill_gradientn(colours=c("#FDE725FF","#22A884FF","#440154FF")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = 'black'),
        legend.key.height = unit(0.4, "cm"),
        legend.position = "top",
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_blank()); p1

# Heatmap for chemotype BThu
p2 <- ggplot(features_subset %>%
               dplyr::filter(Chemotype == "BThu") %>% 
               mutate(Treatment = forcats::fct_recode(Treatment,
                                                      "Mock-treated" = "Mock_treated",
                                                      "Low-dose" = "Low_dose",
                                                      "High-dose" = "High_dose")),
             aes(y = bucket_id, x = Treatment, fill = Intensity_scaled)) + 
  geom_tile() +
  labs(y = NULL,
       fill = "Standardized and scaled feature intensity") +
  facet_grid(Fraction ~ Chemotype, scales = "free") +
  scale_fill_gradientn(colours=c("#FDE725FF","#22A884FF","#440154FF")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.key.height = unit(0.4, "cm"),
        legend.position = "top",
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)); p2

# Combine heatmaps with common legend
ggarrange(p1, p2,
          ncol=2, nrow=1,
          labels = c("(A)", "(B)"),
          common.legend = TRUE,
          legend="top") +
  theme_cowplot()

# Save
if(!dir.exists(here("Figures"))){
  dir.create(here("Figures"))
}
ggsave(plot = last_plot(),
       filename = here("Figures", "Figure_2.tiff"),
       dpi = 600,
       width = 20,
       height = 24,
       units = "cm",
       bg = "transparent")

# =========================T=H=E==E=N=D====================================== #
