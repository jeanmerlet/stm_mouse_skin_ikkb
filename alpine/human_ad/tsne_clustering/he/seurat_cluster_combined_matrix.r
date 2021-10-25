library(Seurat)

mtx_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/combined_human-ad_mtx.tsv'

mtx <- read.table(file=mtx_path, header=TRUE, row.names=1, sep='\t')

# convert to seurat object
mtx <- CreateSeuratObject(counts=mtx)
all_genes <- rownames(mtx)

# preprocessing stuff
mtx <- FindVariableFeatures(mtx, nfeatures=2000)
mtx <- ScaleData(mtx, features=all_genes)
mtx <- RunPCA(mtx, features=VariableFeatures(object=mtx))
mtx <- JackStraw(mtx, num.replicate=100) # this step takes ~15 min
mtx <- ScoreJackStraw(mtx, dims=1:20)

# clustering
#mtx <- FindNeighbors(mtx, dims=1:20)
mtx <- FindNeighbors(mtx, dims=1:10)
mtx <- FindClusters(mtx, resolution=0.5)

out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/clusters_10dim.tsv'
write.table(as.data.frame(Idents(mtx)), file=out_path, sep='\t', quote=F)
