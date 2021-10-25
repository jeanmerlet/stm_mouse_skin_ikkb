suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scran)
  library(cellassign)
  #library(pheatmap)
})
reticulate::use_condaenv("/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rhea/envs/cellassign")

data_dir <- "/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/"

ct_suffix <- "_cell_types_mod1.tsv"
params_suffix <- "_cell_mle_params_mod1.Rda"

umi_norm_by_cell <- function(exp_mat) {
  norm_factors <- colSums(exp_mat)
  exp_mat <- t(t(exp_mat) / norm_factors * 10^4)
  return(exp_mat)
}

experiments <- list('FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848')
for (experiment in experiments) {
celltypes_dir <- paste0(data_dir, experiment, "/cell_types/")
name <- substr(experiment, 11, 16)
ct_path <- paste0(celltypes_dir, name, ct_suffix)
celltypes <- read.table(ct_path, sep="\t", stringsAsFactors=FALSE, colClasses=c('character', 'character'))

matrix_dir <- paste0(data_dir, experiment, "/filtered_feature_bc_matrix")
sse <- read10xCounts(matrix_dir, col.names=T, type="sparse", version="3")
exp_mat <- assays(sse)$counts
exp_mat <- umi_norm_by_cell(exp_mat)
barcodes <- paste(colnames(sse), name, sep="_")
gene_symbols <- toupper(rowData(sse)$Symbol)
rownames(exp_mat) <- gene_symbols

t_bc <- celltypes[celltypes[, 2] == 'T cells', 1]
t_bc_idx <- barcodes %in% t_bc
t_genes <- c("CD3G", "TRBC1", "TRAC", "CD4", "CD8A", "TRDV4")
t_gene_idx <- match(t_genes, gene_symbols)
t_exp_mat <- exp_mat[t_gene_idx, t_bc_idx] 
write.table(as.matrix(t_exp_mat), file=paste0("t_cells_", name, ".tsv"), row.names=TRUE, col.names=NA, sep="\t")
i <- 0
for (r in 1:nrow(t_exp_mat)) {
  i <- i + 1
  print(t_genes[i])
  row <- t_exp_mat[r, ]
  print(sum(row))
  #print(max(row))
  #print(mean(row))
  #print(median(row))
}
}
break

ilc1_bc <- celltypes[celltypes[, 2] == 'ILC1', 1]
ilc1_bc_idx <- barcodes %in% ilc1_bc
ilc1_genes <- c("KLRD1", "IL2RB", "KLRK1", "NKG7", "CTSW", "CCL5", "SERPINB9B", "NCR1", "KLRB1", "TBX21", "CD3G", "CD4", "PTPRC")
ilc1_gene_idx <- match(ilc1_genes, gene_symbols)
ilc1_exp_mat <- exp_mat[ilc1_gene_idx, ilc1_bc_idx] 
write.table(as.matrix(ilc1_exp_mat), file=paste0("ilc1_cells_", name, ".tsv"), row.names=TRUE, col.names=NA, sep="\t")
i <- 0
for (r in 1:nrow(ilc1_exp_mat)) {
  i <- i + 1
  print(ilc1_genes[i])
  row <- ilc1_exp_mat[r, ]
  print(sum(row))
  #print(max(row))
  #print(mean(row))
  #print(median(row))
}

ilc2_bc <- celltypes[celltypes[, 2] == 'ILC2', 1]
ilc2_bc_idx <- barcodes %in% ilc2_bc
ilc2_genes <- c("IL7R", "GATA3", "IL1RL1", "LY6A", "ICOS", "CD3G", "CD4", "PTPRC")
ilc2_gene_idx <- match(ilc2_genes, gene_symbols)
ilc2_exp_mat <- exp_mat[ilc2_gene_idx, ilc2_bc_idx] 
write.table(as.matrix(ilc2_exp_mat), file=paste0("ilc2_cells_", name, ".tsv"), row.names=TRUE, col.names=NA, sep="\t")
i <- 0
for (r in 1:nrow(ilc2_exp_mat)) {
  i <- i + 1
  print(ilc2_genes[i])
  row <- ilc2_exp_mat[r, ]
  print(sum(row))
  #print(max(row))
  #print(mean(row))
  #print(median(row))
}

gbt_bc <- celltypes[celltypes[, 2] == 'Gamma delta T cells', 1]
gbt_bc_idx <- barcodes %in% gbt_bc
gbt_genes <- c("NKG7", "KLRD1", "GZMA", "CCL5", "TRDC", "CD8A", "PTPRC")
gbt_gene_idx <- match(gbt_genes, gene_symbols)
gbt_exp_mat <- exp_mat[gbt_gene_idx, gbt_bc_idx] 
write.table(as.matrix(gbt_exp_mat), file=paste0("gbt_cells_", name, ".tsv"), row.names=TRUE, col.names=NA, sep="\t")
i <- 0
for (r in 1:nrow(gbt_exp_mat)) {
  i <- i + 1
  print(gbt_genes[i])
  row <- gbt_exp_mat[r, ]
  print(sum(row))
  #print(max(row))
  #print(mean(row))
  #print(median(row))
}

dc_bc <- celltypes[celltypes[, 2] == 'Dendritic cells', 1]
dc_bc_idx <- barcodes %in% dc_bc
dc_genes <- c("BATF3", "CXCL16", "CTSS", "IRF8", "H2-AB1", "H2-EB1", "ITGAX", "CD3G", "CD19", "LY6G", "LY6C1", "LY6C2", "PTPRC")
dc_gene_idx <- match(dc_genes, gene_symbols)
dc_exp_mat <- exp_mat[dc_gene_idx, dc_bc_idx] 
write.table(as.matrix(dc_exp_mat), file=paste0("dc_cells_", name, ".tsv"), row.names=TRUE, col.names=NA, sep="\t")
i <- 0
for (r in 1:nrow(dc_exp_mat)) {
  i <- i + 1
  print(dc_genes[i])
  row <- dc_exp_mat[r, ]
  print(sum(row))
  #print(max(row))
  #print(mean(row))
  #print(median(row))
}
}
