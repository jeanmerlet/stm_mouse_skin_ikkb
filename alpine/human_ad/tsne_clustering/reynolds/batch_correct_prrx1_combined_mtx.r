library(Seurat)

mtx_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/combined_all_prrx1_pos.tsv'

mtx <- read.table(file=mtx_path, header=TRUE, row.names=1, sep='\t')
mtx <- CreateSeuratObject(counts=mtx)

# create case-control id list
conditions <- c(rep('control', 12360), rep('ad', 14296))
names(conditions) <- colnames(mtx)
mtx <- AddMetaData(object=mtx, metadata=conditions, col.name='condition')
condition_list <- SplitObject(mtx, split.by='condition')

for (i in 1:length(condition_list)) {
    condition_list[[i]] <- FindVariableFeatures(condition_list[[i]], selection.method='vst', nfeatures=2000, verbose=FALSE)
}

# batch correction
reference_list <- condition_list[c('control', 'ad')]
mtx_anchors <- FindIntegrationAnchors(object.list=reference_list, dims=1:30)
mtx <- IntegrateData(anchorset=mtx_anchors, dims=1:30)

# write batch-corrected matrix to file
out_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/batch_combined_all_prrx1_pos.tsv'
write.table(mtx@assays$integrated@data, out_path, sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
