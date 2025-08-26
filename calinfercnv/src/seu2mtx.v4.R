# 2025-08-26 v4.2

library(Seurat)
library(dplyr)
library(yaml)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
print(args)

dyml <- yaml.load_file(args[1])
rds <- dyml$rds
assay <- dyml$assay
ctl <- dyml$ctl
omit <- dyml$omit

ref_gene <- dyml$ref_gene
group <- dyml$group
samples <- dyml$samples
outdir <- dyml$outdir

dir.create(outdir)


seurat_obj <- readRDS(rds)
metadata <- seurat_obj@meta.data



cell_anno <- metadata[,c('seurat_clusters', 'orig.ident')]
colnames(cell_anno) <- c('seurat_clusters', 'sample')
cell_anno$cell <- rownames(metadata)
if (nchar(omit) > 0){
    omit_group <- strsplit(omit, ',')[[1]]
    print(omit_group)
    cell_anno <- cell_anno %>% filter(! seurat_clusters %in% omit_group)
    print(unique(cell_anno$seurat_clusters))
}


cell_anno$group <- ifelse(group=='sample', cell_anno$sample, cell_anno$seurat_clusters)
cell_anno$group <- as.character(cell_anno$group)
ref_group <- strsplit(ctl, ',')[[1]]
print(ref_group)
cell_anno$group[cell_anno$seurat_clusters %in% ref_group] <- 'Control'


if (samples != 'all'){
    samples <- stringr::str_split(samples, ',')[[1]]
    cell_anno <- cell_anno %>% filter(cell_anno$group %in% c('Control', samples))
}


fwrite(cell_anno, file=paste0(outdir, '/metadata.used.tsv'), sep='\t', row.names=F)


cell_anno <- cell_anno[, c('cell', 'group')]
fwrite(cell_anno, file=paste0(outdir, '/cellInfo.txt'), sep='\t', row.names=F, col.names=F)


data_genes <- row.names(seurat_obj@assays[[assay]]$counts)

d_gene <- read.table(ref_gene, sep='\t', header=F)
d_gene <- d_gene %>% distinct(V1, .keep_all = TRUE)
rownames(d_gene) <- d_gene$V1
d_gene$V1 <- NULL
used_genes <- intersect(data_genes, rownames(d_gene))
d_gene2 <- d_gene[used_genes,]
fwrite(d_gene2, file=paste0(outdir, '/geneLoci.txt'), sep = "\t", quote=F, row.names=T, col.names=F)


assay_count <- seurat_obj@assays[[assay]]$counts[used_genes, cell_anno$cell]
count_mtx <- as.data.frame(as.matrix(assay_count))
rownames(count_mtx) <- used_genes
fwrite(count_mtx, file=paste0(outdir, '/expCount.txt'), sep = "\t", quote=F, row.names=T, col.names=T)





