# INFO ####
# This R script sources R_Figure_1a_and_1b.R, R_Figure_1c_and_1d.R and
# R_Figure_1e_and_1f.R to create Figure 1
# Author: Silvia Eckert
# Generated: 2024-04-25

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "microdiluteR", "ggplot2", "dplyr"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Analysis_BMA/R_Figure_1.R") # Set root directory
setwd(here::here())
getwd()
library(microdiluteR)
library(ggplot2)
library(dplyr)

# =========================================================================== #

# SOURCE R SCRIPTS ####
# R_Figure_1a_and_1b.R
source(here::here("Analysis_BMA", "R_Figure_1a_and_1b.R"))
# R_Figure_1a_and_1b.R
source(here::here("Analysis_BMA", "R_Figure_1c_and_1d.R"))
# R_Figure_1a_and_1b.R
source(here::here("Analysis_BMA", "R_Figure_1e_and_1f.R"))

# =========================================================================== #

# GENERATE FIGURE ####
ggpubr::ggarrange(Figure_1a, Figure_1b,
                  Figure_1c, Figure_1d,
                  Figure_1e, Figure_1f,
                  ncol = 2,
                  nrow = 3,
                  heights = c(0.85, 0.85, 1),
                  labels = paste0("(", LETTERS[1:6], ")"),
                  common.legend = TRUE,
                  legend = "top")

# Save
if(!dir.exists(here("Figures"))){
  dir.create(here("Figures"))
}
ggsave(plot = last_plot(),
       filename = here("Figures", "Figure_1.tiff"),
       dpi = 600,
       width = 20,
       height = 24,
       units = "cm",
       bg = "transparent")

# =========================T=H=E==E=N=D====================================== #
