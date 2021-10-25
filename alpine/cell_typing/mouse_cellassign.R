# needed for cellassign to correctly locate tensorflow
reticulate::import("tensorflow")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scran)
  library(cellassign)
})

# create a named list with all genes associated with each cell type
markers <- list()

# Fibroblasts #
markers[["Dermal fibroblasts 1"]] <- c("SPARC", "COL1A1")
markers[["Dermal fibroblasts 2"]] <- c("DCN", "LUM")
markers[["Hypodermal fibroblasts 1"]] <- c("CYGB", "GPX3")
markers[["Hypodermal fibroblasts 2"]] <- c("ANXA3", "PLAC8")
markers[["Dermal sheath"]] <- c("LRRC15", "ENPP2")
markers[["Dermal papilla"]] <- c("CRABP1", "RASD1")
markers[["Myofibroblasts"]] <- c("MYL9", "ACTA2", "AOC3")

# Other
markers[["Keratinocytes"]] <- c("KRT5", "KRT14", "DSP")
markers[["Smooth muscle cells"]] <- c("ACTA2", "TPM2", "MYH11", "MYLK")
markers[["Skeletal muscle"]] <- c("DES", "ACTA1")
markers[["Schwann cells"]] <- c("MBP", "CNP")
markers[["Endothelial cells"]] <- c("GPIHBP1", "FABP4", "LY6C1")

# Macrophages
markers[["Mono / Macro"]] <- c("CD14", "ADGRE1")

# Lymphocytes
markers[["T cells"]] <- c("CD3G", "NKG7", "THY1")
markers[["T reg"]] <- c("FOXP3", "CTLA4", "TNFRSF4")
markers[["gd T cels"]] <- c("TRDC", "TRDV4", "CXCR6")
markers[["NKT cells"]] <- c("CD3G", "NKG7", "CD8A")
markers[["B cells"]] <- c("CD19", "MS4A1", "CD79A")
markers[["Plasma cells"]] <- c("IGHM", "IGKC")
markers[["NK cells"]] <- c("NKG7", "FCER1G", "KLRB1")
markers[["ILC2"]] <- c("AREG", "GATA3", "IL13", "IL5")

# Other immune cells
markers[["Langerhans cells"]] <- c("CD207", "MFGE8")
markers[["Dendritic cells"]] <- c("CST3", "CLEC9A", "CD207")
markers[["Mast cells"]] <- c("CPA3", "CMA1", "MCPT4", "TPSB2")
markers[["Granulocytes"]] <- c("CSF3", "CSF3R", "IL1R2")

# save celltype-marker list for future reference
marker_meta_dir <- "/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rna-seq_tools/cell_typing/cellassign/rhea/marker_gene_lists"
marker_meta_path <- paste0(marker_meta_dir, "/marker_genes_ct_jean-kang1.txt")
capture.output(markers, file=marker_meta_path)

# convert into binary gene by celltype matrix
marker_gene_matrix <- marker_list_to_mat(markers)
print('done setting up marker genes')

experiments <- list('FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848')
experiments <- list('FGC2063_5_86846', 'FGC2063_5_86847')
for (experiment in experiments) {
matrix_dir <- paste0("/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/", experiment, "/filtered_feature_bc_matrix")
### gene count matrix ###
gene_count_matrix <- read10xCounts(matrix_dir, col.names = TRUE, type = "sparse", version = "3")
# compute size factors (normalize) using ***all*** genes
gene_count_matrix <- computeSumFactors(gene_count_matrix)
# change from ENS to SYMBOL to match marker genes
rownames(gene_count_matrix) <- toupper(rowData(gene_count_matrix)$Symbol)
# remove genes not expressed by any cells
gene_count_matrix <- gene_count_matrix[rowSums(assays(gene_count_matrix)$counts) != 0, ]

out_dir <- paste0("/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/", experiment, "/cell_types/")
intersection <- intersect(rownames(gene_count_matrix), rownames(marker_gene_matrix))
current_matrix <- gene_count_matrix[intersection, ]
print(setdiff(rownames(current_matrix), rownames(marker_gene_matrix)))
current_genes <- marker_gene_matrix[intersection, ]
# remove cells with no mapping counts to the current marker genes in the current gene count matrix
current_matrix <- current_matrix[, colSums(assays(current_matrix)$counts) != 0]
### cellassign fit ###
print('starting fit...')
start_time <- proc.time()
fit <- cellassign(exprs_obj = current_matrix,
                  marker_gene_info = current_genes,
                  s = sizeFactors(current_matrix),
                  learning_rate = 0.01,
                  verbose = TRUE)
end_time <- proc.time() - start_time
print('done fitting')
print(end_time)
print(fit)

name <- substr(experiment, 11, 16)
types_path <- paste0(out_dir, name, "_cell_types_jean-kang1.tsv")
mle_path <- paste0(out_dir, name, "_cell_mle_params_jean-kang1.Rda")

cell_type <- fit$cell_type
mle_params <- fit$mle_params
barcodes <- paste(colnames(current_matrix), name, sep="_")
names(cell_type) <- barcodes

write.table(cell_type, types_path, sep="\t", col.names=FALSE)
save(mle_params, file=mle_path)
}
