suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scran)
  library(cellassign)
  library(pheatmap)
})
reticulate::use_condaenv("/Users/6j9/anaconda3/envs/cellassign")

data_dir <- "/Users/6j9/projects/mouse/cytoscape/coloring/raw/"
mcl_dir <- "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/"
out_dir <- "/Users/6j9/projects/mouse/plots/out/celltype_heatmaps/"
fname <- paste0(out_dir, "kang_all.pdf")

mcl_path <- paste0(mcl_dir, "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.20.txt")
mcl_clusters <- strsplit(readLines(mcl_path), split="\t")

barcodes_name <- "_cell_types_kang-split.tsv"
params_name <- "_cell_mle_params_kang-split.Rda"
cellprobs <- data.frame()
names <- c()
experiments <- list('FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848')
for (experiment in experiments) {
  name <- substr(experiment, 11, 16)
  bc_path <- paste0(data_dir, name, barcodes_name)
  mle_path <- paste0(data_dir, name, params_name)
  barcodes <- read.table(bc_path, sep="\t", stringsAsFactors=FALSE, colClasses=c('character', 'NULL'))
  load(mle_path)
  cellprobs <- rbind(cellprobs, data.frame(mle_params$gamma, row.names=unlist(barcodes)))
}

#clusters <- c(1)
sub_cellprobs <- data.frame()
for (cluster in clusters){
  sub_cellprobs <- rbind(sub_cellprobs, cellprobs[mcl_clusters[[cluster]], ])
}
colnames(sub_cellprobs)[colnames(sub_cellprobs) == "Erythroid.like.and.erythroid.precursor.cells"] <- "Erythroid-like"
colnames(sub_cellprobs)[colnames(sub_cellprobs) == "Plasmacytoid.dendritic.cells"] <- "Plasmacytoid"
colMax <- function(data) sapply(data, max, na.rm = TRUE)
#sub_cellprobs <- sub_cellprobs[ ,colMax(sub_cellprobs) > 0.5]
#sub_cellprobs <- sub_cellprobs[ ,colSums(sub_cellprobs) >= 10]
colnames(sub_cellprobs) <- gsub("[.]", " ", colnames(sub_cellprobs))
colnames(cellprobs) <- gsub("[.]", " ", colnames(sub_cellprobs))
colnames(cellprobs)[colnames(sub_cellprobs) == "Prrx1  fibroblasts"] <- "Prrx1+"
colnames(cellprobs)[colnames(sub_cellprobs) == "fibrogenic papillary fibroblasts"] <- "Fibrogenic"
colnames(cellprobs)[colnames(sub_cellprobs) == "hypodermal fibroblasts"] <- "Hypodermal"
colnames(cellprobs)[colnames(sub_cellprobs) == "reticular fibroblasts"] <- "Reticular"
colnames(cellprobs)[colnames(sub_cellprobs) == "monocyte macrophage chemokine"] <- "MM chemokine"
colnames(cellprobs)[colnames(sub_cellprobs) == "eosinophil chemokine"] <- "E chemokine"

#pheatmap(sub_cellprobs, filename=fname, show_rownames=FALSE, treeheight_row=0, treeheight_col=0, fontsize=16, fontsize_col=10, width=15, height=5, angle=90)
pheatmap(cellprobs, filename=fname, show_rownames=FALSE, treeheight_row=0, treeheight_col=0, fontsize=16, fontsize_col=10, width=18, height=5, angle=90)
