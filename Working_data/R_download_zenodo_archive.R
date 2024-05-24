# INFO ####
# This R script downloads the raw data underlying the corresponding
# manuscript from the zenodo repository
# Author: Silvia Eckert
# Generated: 2024-05-13

# =========================================================================== #

# PACKAGES ####
# Uncomment to install packages if necessary
# install. packages(c("here", "zen4R"))
# Set working directory with .here file in root directory of project
library(here)
here::i_am("Working_data/R_download_zenodo_archive.R") # Set root directory
setwd(here::here())
getwd()
library(zen4R)

# =========================================================================== #

# FETCH DATA ####
# Download archive
doi <- ""
download_zenodo(doi,
                path = here("Working_data"),
                logger = "INFO")
# Unzip archive
zip_file_path <- list.files(here("Working_data"),
                            full.names = T,
                            pattern = ".zip"); zip_file_path
for (path in zip_file_path) {
  utils::unzip(path,
               exdir = here("Working_data"))
}

# UNZIP FOLDERS ####
zip_folder_path <- list.files(here("Working_data"),
                              full.names = T,
                              pattern = "plates.zip"); zip_folder_path
for (folder in zip_folder_path) {
  utils::unzip(folder,
               exdir = here("Working_data"))
}

# =========================T=H=E==E=N=D====================================== #
