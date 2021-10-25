library(Seurat)

mtx_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/combined_human-ad_mtx.tsv'
mtx <- read.table(file=mtx_path, header=TRUE, row.names=1, sep='\t')

# convert to seurat object
mtx <- CreateSeuratObject(counts=mtx)

# create case-control id list
conditions <- c(rep('control', 19537), rep('ad', 33452))
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
out_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/batch-corrected_combined_human-ad_mtx.tsv'
write.table(mtx@assays$integrated@data, out_path, sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)

# rest of preprocessing stuff
mtx <- read.table(file=out_path, header=TRUE, row.names=1, sep='\t')
mtx <- CreateSeuratObject(counts=mtx)
all_genes <- rownames(mtx)

mtx <- FindVariableFeatures(mtx, nfeatures=2000) # need to run this for PCAs to work
mtx <- ScaleData(mtx, features=all_genes)
mtx <- RunPCA(mtx, features=VariableFeatures(object=mtx))
mtx <- JackStraw(mtx, num.replicate=100) # this step takes ~8 min
mtx <- ScoreJackStraw(mtx, dims=1:20)

# clustering
mtx <- FindNeighbors(mtx, dims=1:20)
#mtx <- FindNeighbors(mtx, dims=1:10)
mtx <- FindClusters(mtx, resolution=0.5)

out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/reynolds/clusters_20dim.tsv'
write.table(as.data.frame(Idents(mtx)), file=out_path, sep='\t', quote=F)
