suppressPackageStartupMessages({
    library(DESeq2)
    library(DropletUtils)
    library(readr)
    library(scran)
    library(apeglm)
    library(stringr)
})

matrix_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/r_combined_all_prrx1_pos_unprocessed.tsv'
out_dir <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/de_results/'
genes_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/all_genes.txt'

# create sce
con <- file(matrix_path, 'r')
all_bcs <- unlist(str_split(readLines(con, n=1), '\t'))
close(con)
counts <- as.matrix(read.table(matrix_path, sep='\t', header=FALSE, skip=1))

# gene names
con <- file(genes_path, 'r')
genes <- readLines(con)
close(con)

get_condition_from_bc <- function(barcode) {
    if (str_detect(barcode, '809') == TRUE) {
        return('ad')
    } else {
        return('control')
    }
}

get_batch_from_bc <- function(barcode) {
    if (str_detect(barcode, '4820') == TRUE) {
        return('batch1')
    } else if (str_detect(barcode, '810') == TRUE) {
        return('batch2')
    } else {
        return('batch3')
    }
}

# create DE design columns
conditions <- character()
batches <- character()
for (bc in all_bcs) {
    conditions <- c(conditions, get_condition_from_bc(bc))
    batches <- c(batches, get_batch_from_bc(bc))
}

meta_data <- data.frame('condition'=conditions, row.names=all_bcs, stringsAsFactors=TRUE)

# create DESeq2 data set
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta_data, design= ~ condition)
rownames(dds) <- genes

# compute scran normalization
scran_clusters <- quickCluster(dds)
dds <- computeSumFactors(dds, clusters=scran_clusters)

# Healthy vs. AD
dds$condition <- factor(dds$condition, levels=c('ad', 'control'))

start_time <- proc.time()
dds <- DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, parallel=FALSE)
print(proc.time() - start_time)

res <- results(dds, alpha=0.05, parallel=TRUE)
writeLines(capture.output(summary(res)), paste0(out_dir, 'summary_healthy_vs_ad.txt'))
ordered_res <- res[order(res$padj), ]
write.table(as.data.frame(ordered_res), sep='\t', file=paste0(out_dir, 'table_healthy_vs_ad.tsv'))
