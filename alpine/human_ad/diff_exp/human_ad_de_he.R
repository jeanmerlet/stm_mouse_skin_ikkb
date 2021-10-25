suppressPackageStartupMessages({
    library(DESeq2)
    library(DropletUtils)
    library(readr)
    library(scran)
    library(apeglm)
    library(stringr)
})

matrix_dir_prefix <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/aligned'
clusters_path <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/clusters_20dim.tsv'
out_dir <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/de_results/'

# create sce
sample_dirs <- list.dirs(path=matrix_dir_prefix, full.names = TRUE, recursive=FALSE)
sample_names <- str_extract(sample_dirs, 'SRR\\d+') 

combine_samples <- function(sample_names) {
    matrix_dirs <- character()
    barcodes <- character()
    for (sample_name in sample_names) {
        matrix_dir <- paste0(matrix_dir_prefix, '/', paste0(sample_name, '_1._Solo.out/Gene/filtered'))
        matrix_dirs <- c(matrix_dirs, matrix_dir)
        bc_path <- paste0(matrix_dir, '/barcodes.tsv')
        bc <- readLines(bc_path)
        barcodes <- c(barcodes, paste(bc, sample_name, sep='_'))
    }
    sce <- read10xCounts(matrix_dirs, col.names=FALSE, type='sparse', version='3', compressed=FALSE)
    colnames(sce) <- barcodes
    rownames(sce) <- toupper(rowData(sce)$Symbol)
    return(sce)
}
sce <- combine_samples(sample_names)
counts <- as.matrix(assay(sce))

bc_to_condition_map <- list('SRR11396159'='AD_LS',
                            'SRR11396160'='AD_LS',
                            'SRR11396161'='AD_NL',
                            'SRR11396162'='Healthy',
                            'SRR11396163'='AD_LS',
                            'SRR11396164'='Healthy',
                            'SRR11396165'='AD_LS',
                            'SRR11396166'='Healthy',
                            'SRR11396167'='Healthy',
                            'SRR11396168'='Healthy',
                            'SRR11396169'='AD_NL',
                            'SRR11396170'='Healthy',
                            'SRR11396171'='Healthy',
                            'SRR11396172'='AD_NL',
                            'SRR11396173'='AD_NL',
                            'SRR11396174'='AD_NL',
                            'SRR11396175'='Healthy')

get_condition_from_bc <- function(barcode) {
    id <- substr(barcode, 18, 29)
    return(bc_to_condition_map[[id]])
}

bc_to_sample_map <- list('S1' = 'SRR11396159',
                         'S2' = 'SRR11396160',
                         'S3' = 'SRR11396161',
                         'S4' = 'SRR11396162',
                         'S5' = 'SRR11396163',
                         'S6' = 'SRR11396164',
                         'S7' = 'SRR11396165',
                         'S8' = 'SRR11396166',
                         'S9' = 'SRR11396167',
                         'S10' = 'SRR11396168',
                         'S11' = 'SRR11396169',
                         'S12' = 'SRR11396170',
                         'S13' = 'SRR11396171',
                         'S14' = 'SRR11396172',
                         'S15' = 'SRR11396173',
                         'S16' = 'SRR11396174',
                         'S17' = 'SRR11396175')

# restrict to mesenchymal cells by KNN cluster
clusters <- read.table(clusters_path, sep='\t', header=FALSE)
colnames(clusters) <- c('bc', 'cluster')
mc_bcs <- clusters[clusters[ ,'cluster'] %in% c(3, 4, 5, 8), ]['bc'][[1]]
i <- 0
for (bc in mc_bcs) {
    barcode <- str_extract(bc, '[A-Z]+$')
    sample <- str_extract(bc, 'S\\d+')
    i <- i + 1
    mc_bcs[[i]] <- paste0(barcode, '_', bc_to_sample_map[[sample]])
}
all_bcs <- colnames(sce)

# restrict to PRRX1+ Mesenchymal cells only
prrx1_counts <- counts['PRRX1', ]
bcs_to_keep <- character()
for (bc in mc_bcs) {
    if (bc %in% names(prrx1_counts)) {
        if (prrx1_counts[bc] > 0) {
            bcs_to_keep <- c(bcs_to_keep, bc)
        }
    }
}

de_idx <- which(all_bcs %in% mc_bcs)
all_bcs <- all_bcs[de_idx]
de_counts <- counts[, de_idx]

# create DE design columns
conditions <- character()
for (bc in all_bcs) {
    conditions <- c(conditions, get_condition_from_bc(bc))
}
meta_data <- data.frame('condition'=conditions, row.names=all_bcs, stringsAsFactors=TRUE)

# create DESeq2 data set
dds <- DESeqDataSetFromMatrix(countData=de_counts, colData=meta_data, design= ~ condition)

# remove lowly expressed genes
keep_idx <- rowSums(counts(dds)) >= 10
dds <- dds[keep_idx, ]

# compute scran normalization
scran_clusters <- quickCluster(dds)
dds <- computeSumFactors(dds, clusters=scran_clusters)

# LS vs. H
dds$condition <- factor(dds$condition, levels=c('AD_LS', 'Healthy'))
dds$condition <- droplevels(dds$condition)
dds <- dds[, dds$condition %in% c("AD_LS", "Healthy") ]

start_time <- proc.time()
dds <- DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, parallel=FALSE)
print(proc.time() - start_time)

res <- results(dds, alpha=0.05, parallel=TRUE)
writeLines(capture.output(summary(res)), paste0(out_dir, 'summary_ls-ad_vs_h.txt'))
ordered_res <- res[order(res$padj), ]
write.table(as.data.frame(ordered_res), sep='\t', file=paste0(out_dir, 'table_ls-ad_vs_h.tsv'))

# LS vs. NL
dds$condition <- factor(dds$condition, levels=c('AD_LS', 'AD_NL'))
dds$condition <- droplevels(dds$condition)
dds <- dds[, dds$condition %in% c("AD_LS", "AD_NL") ]

start_time <- proc.time()
dds <- DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, parallel=FALSE)
print(proc.time() - start_time)

res <- results(dds, alpha=0.05, parallel=TRUE)
writeLines(capture.output(summary(res)), paste0(out_dir, 'summary_ls-ad_vs_nl.txt'))
ordered_res <- res[order(res$padj), ]
write.table(as.data.frame(ordered_res), sep='\t', file=paste0(out_dir, 'table_ls-ad_vs_nl.tsv'))
