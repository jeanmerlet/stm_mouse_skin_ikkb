suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scran)
  library(cellassign)
})
print('done loading libraries')
reticulate::use_condaenv("cellassign")

matrix_dir = "/Users/6j9/projects/mouse/data/FGC2091_7_92848/filtered_feature_bc_matrix/"
markers_filepath = "/Users/6j9/projects/mouse/annotation/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv"

#experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']

### marker gene matrix ###
markers <- read.table(markers_filepath, h=T, sep='\t', quote='', stringsAsFactors=FALSE)
# remove non-canonical markers
markers <- markers[!is.na(markers$canonical.marker), ]
# remove human-only markers
markers <- markers[markers$species != "Hs", ]
# create a named list with all genes associated with each cell type
marker_gene_list <- list()
for (cell_type in unique(markers$cell.type)) {
  marker_gene_list[[cell_type]] <- markers[markers$cell.type == cell_type, ][ ,"official.gene.symbol"]
}
# convert into binary gene by celltype matrix
marker_gene_matrix <- marker_list_to_mat(marker_gene_list)
print('done setting up marker genes')

### gene count matrix ###
gene_count_matrix = read10xCounts(matrix_dir, col.names = TRUE, type = "sparse", version = "3")
# compute size factors using ***all*** genes
gene_count_matrix <- computeSumFactors(gene_count_matrix)
# change from ENS to SYMBOL to match marker genes
rownames(gene_count_matrix) <- toupper(rowData(gene_count_matrix)$Symbol)
print('done setting up gene count matrix')

# subset to ***only marker genes*** for use by cellassign
intersection <- intersect(rownames(gene_count_matrix), rownames(marker_gene_matrix))
gene_count_matrix <- gene_count_matrix[intersection, ]
marker_gene_matrix <- marker_gene_matrix[intersection, ]
# remove cells with no mapping counts
expression_matrix <- assays(gene_count_matrix)$counts
gene_count_matrix <- gene_count_matrix[ ,colSums(expression_matrix) != 0]
# remove genes with no mapping counts
gene_count_matrix <- gene_count_matrix[rowSums(expression_matrix) != 0, ]
marker_gene_matrix <- marker_gene_matrix[rownames(gene_count_matrix), ]

print(dim(gene_count_matrix))
print(dim(marker_gene_matrix))
break

# assign correct dimension size factors (computed earlier) for use by cellassign
s <- sizeFactors(gene_count_matrix)

fit <- cellassign(exprs_obj = gene_count_matrix,
                  marker_gene_info = marker_gene_matrix,
                  s = s,
                  learning_rate = 1e-2,
                  verbose = TRUE)

#save celltypes(fit)
#save cellprobs(fit)

print(fit)
