

####----DRPPM - PATH - SURVEIOR R Package Installation----####


##--This Will install immunedeconv package--##
##---R Version 4.1 or greater is required---##


## Check if packages are installed
packages <- c("shiny","shinythemes","shinyjqui","gtsummary","tidyr","RColorBrewer","immunedeconv",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap","stringr",
              "plotly","readr","shinycssloaders","survminer","gridExtra","viridis")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

## Check if Bioconductor specific packages are installed
bioCpacks <- c("GSVA","clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))
