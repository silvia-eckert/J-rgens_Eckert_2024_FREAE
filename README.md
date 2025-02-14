
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R code from: Jürgens and Eckert (2024)

![](https://img.shields.io/badge/repo%20status-archived-orange.svg)
[![DOI](https://zenodo.org/badge/805307619.svg)](https://zenodo.org/doi/10.5281/zenodo.12825708)

# :open_file_folder: Repository Structure

- **Working_data/:** Contains the initial R script to download and unzip
  the data archived at [zenodo](https://zenodo.org/)
- **Analysis_BMA/:** Contains all R scripts to analyse the broth
  microdilution assay (BMA) dataset
- **Analysis_LCMS/:** Contains all R scripts to analyse the liquid
  chromatography-mass spectrometry (LCMS) dataset

# :package: Dependencies

To be able to smoothly reproduce the figures and tables, make sure you
have [R v4.3.3](https://cran.uni-muenster.de/index.html) and the
following R packages installed:

``` r
install.packages(c("here", "microdiluteR", "dplyr", "tidyr", "ggplot2", "ggpubr", "cowplot", "vegan", "reshape2"))
```

# :arrow_forward: Quick start

To reproduce the analyses and visualizations from the manuscript in a
structured and guided manner:

1.  **Download** [R v4.3.3](https://cran.uni-muenster.de/index.html) and
    [RStudio](https://posit.co/download/rstudio-desktop/).
2.  **Clone the repository** via
    [Git](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository),
    [VS
    code](https://learn.microsoft.com/en-gb/azure/developer/javascript/how-to/with-visual-studio-code/clone-github-repository?tabs=activity-bar)
    or
    [RStudio](https://argoshare.is.ed.ac.uk/healthyr_book/clone-an-existing-github-project-to-new-rstudio-project.html).
3.  **Keep the structure** as it is, don’t rename the folders or R
    scripts that are already present.
4.  **Open the `.Rproj` file** in RStudio to initiate the cloned
    repository as an RStudio project. This will enable all R scripts to
    automatically detect the root directory of the project.
5.  **Run the `R_download_zenodo_archive.R` script** with `CTRL+ALT+R`
    to download the stored data from the corresponding
    [zenodo](https://zenodo.org/) repository. You don’t have to change
    the settings, since the data will be automatically downloaded and
    unzipped into the `Working_data` folder.
6.  **Analyze BMA data** by running the `R_Figure_1.R` and
    `R_Table_S2.R` scripts. Two new folders `Figures` and `Tables` will
    emerge storing the corresponding figures and tables shown in the
    manuscript.
7.  **Analyze LCMS data** by running the scripts `R_tidy_LCMS_data.R`,
    `"R_Figure_2.R` and `R_Tables_1_and_S3.R` sequentially. Additional
    figures and tables will be added to the folders `Figures` and
    `Tables`.

# :pushpin: Citation

**Citing this repository:**

> Jürgens, J., & Eckert, S. (2024). R code from: Foliar 5-azacytidine
> treatment of Tanacetum vulgare chemotypes alters in vitro growth
> performance of the necrotrophic fungus Botrytis cinerea (v1.0.0)
> \[Computer software\].
> <https://github.com/silvia-eckert/J-rgens_Eckert_2024_FREAE>

**Citing the corresponding manuscript:**

> Jürgens, J., & Eckert, S. (2024). Foliar 5-azacytidine treatment of
> *Tanacetum vulgare* chemotypes alters in vitro growth performance of
> the necrotrophic fungus *Botrytis cinerea*. Frontiers in Epigenetics
> and Epigenomics (submitted).

**Citing the dataset deposited on Zenodo:**

> Jürgens, J., & Eckert, S. (2024). Data from: Foliar 5-azacytidine
> treatment of *Tanacetum vulgare* chemotypes alters in vitro growth
> performance of the necrotrophic fungus *Botrytis cinerea* (v1.0.0)
> \[Data set\]. *Zenodo*.
> [https://doi.org/10.5281/zenodo.11281235](https://www.doi.org/10.5281/zenodo.11281235)

# :scroll: License

This repository and the files therein are licensed under the [GPL-3.0
license](https://www.gnu.org/licenses/gpl-3.0.html.en).
