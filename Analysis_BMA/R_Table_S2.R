# INFO ####
# This R script sources R_Figure_2a_and_2b.R, R_Figure_2c_and_2d.R and
# R_Figure_2e_and_2f.R to create Table 1
# Author: Silvia Eckert
# Generated: 2024-04-25

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "microdiluteR", "ggplot2", "dplyr"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Analysis_BMA/R_Table_S2.R") # Set root directory
setwd(here::here())
getwd()
library(microdiluteR)
library(ggplot2)
library(dplyr)

# =========================================================================== #

# SOURCE R SCRIPTS ####
# R_Figure_2a_and_2b.R
source(here::here("Analysis_BMA", "R_Figure_1a_and_1b.R"))
# R_Figure_2a_and_2b.R
source(here::here("Analysis_BMA", "R_Figure_1c_and_1d.R"))
# R_Figure_2a_and_2b.R
source(here::here("Analysis_BMA", "R_Figure_1e_and_1f.R"))

# =========================================================================== #

# GENERATE SUMMARY ####
# Combined stats in a list
stats_list <- list(Aketo_0uM5AZA_stats = Aketo_0uM5AZA_stats,
                   Aketo_50uM5AZA_stats = Aketo_50uM5AZA_stats,
                   Aketo_100uM5AZA_stats = Aketo_100uM5AZA_stats,
                   BThu_0uM5AZA_stats = BThu_0uM5AZA_stats,
                   BThu_50uM5AZA_stats = BThu_50uM5AZA_stats,
                   BThu_100uM5AZA_stats = BThu_100uM5AZA_stats)
# Subset data and combine to a data frame
table_S2 <- data.frame()
for (stats in names(stats_list)) {
  
  # Extract stats element and subset columns
  stats_subset <- stats_list[[stats]] %>% 
    dplyr::select(c(Group, Experiment, Treatment, Concentration, mean, stderr, n, statistic, df, p, p.signif)) %>% 
    mutate(p.signif = if_else(n < 6, "", p.signif),
           p = if_else(n < 6, NA, p),
           statistic = if_else(n < 6, NA, statistic),
           df = if_else(n < 6, NA, df),
           p = round(p, 3)) %>% 
    dplyr::rename(Chemotype = Group,
                  Fraction = Treatment,
                  Mean = mean,
                  SE = stderr,
                  Conc = Concentration,
                  Treatment = Experiment,
                  Notation = p.signif,
                  S = statistic) %>% 
    mutate(Fraction = factor(Fraction, levels = c("F10", "F18", "F30", "F100"))) %>% 
    dplyr::relocate(Chemotype, Treatment, Fraction, Conc, Mean, SE, n, S, df, p, Notation) %>% 
    dplyr::arrange(Chemotype, Treatment, Fraction)
  
  table_S2 <- rbind(table_S2, stats_subset)
}; table_S2

# save
if(!dir.exists(here("Tables"))){
  dir.create(here("Tables"))
}
write.table(table_S2,
            file = here("Tables",
                        "Table_S2.txt"),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t",
            na = "")

# =========================T=H=E==E=N=D====================================== #
