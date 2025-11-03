install.packages(
c(
  'optparse'
  ),
 Ncpus=6
)

reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat",
	   "cli", "crayon", "R6", "rlang", "SeuratObject", "stringi", "withr")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "ggplot2", "gridExtra", "magrittr", "ggdendro", 
           "here", "rprojroot", "tools", "rjson", "anndata")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("zellkonverter")
BiocManager::install("SingleCellExperiment")
BiocManager::install("rhdf5")
BiocManager::install("scuttle")
BiocManager::install("HDF5Array")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# library(reticulate)
# reticulate::py_install("anndata")

install.packages("devtools", repos = "http://cran.rstudio.com/")

tryCatch({
    devtools::install_github("SGDDNB/ShinyCell")
},
    error = function(e) {
    message("Error installing SGDDNB/ShinyCell")
    message(e$message)
    quit("no", -1)
    }
)

devtools::install_github("cellgeni/sceasy")

#memory.limit(size=2500)