library(ShinyCell)
library(rjson)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(scuttle)
library(zellkonverter)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript make_shinycell.R <adata_object> <ref>\n")
  quit(status = 1)
}
inpFile <- args[1]
tissue <- args[2]
metadata_file <- args[3]

mapping_file <- "/opt/ensembl_to_symbol.json"

options(future.globals.maxSize = 8000 * 1024^2)

sce <- readH5AD(inpFile, use_hdf5 = TRUE)

sce <- logNormCounts(sce, assay.type="unscaled")

#read the mapping file
hugoMapping <- fromJSON(file=mapping_file)

#get the ensemble ids
ensIds <- rownames(sce)

#convert ensemble ids to hugo sybmols
hugos<-unname(hugoMapping[ensIds])

#replace null hugos with the ensId
for (i in 1:length(hugos)){
  h <- unlist(hugos[i])
  if(is.null(h)){
    hugos[i]<-ensIds[i]}
  else{
    hugos[i]<-h
  }
}
#convert to 1-d list
hugos<-unlist(hugos)

#replace the rownames of the sce
rownames(sce)<-hugos

tryCatch({
  scConf = createConfig(sce)
},
error = function(e) {
  reticulate::py_last_error()
  message("Error creating config")
  message(e$message)
  quit("no", -1)
}
)
mainDir = "shinyApps"
subDir = tissue

json_data <- fromJSON(file=metadata_file)
uuid = json_data['Data Product UUID']

tissueDir <- file.path(mainDir, subDir)
shinyDir <- file.path(mainDir, subDir, uuid)
if (!dir.exists(mainDir)) {dir.create(mainDir)}
if (!dir.exists(tissueDir)) {dir.create(tissueDir)}
if (!dir.exists(shinyDir)) {dir.create(shinyDir)}  

options(error = function() {
  sink(stderr())
  on.exit(sink(NULL))
  traceback()
})

title = sprintf("Shiny Cell h5ad % s", tissue)
makeShinyApp(sce, scConf, shiny.dir = shinyDir, shiny.title = title) 
 