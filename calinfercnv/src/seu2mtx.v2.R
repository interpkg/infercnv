# 2025-06-20 v3.2

library(Seurat)
library(dplyr)
library(yaml)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
print(args)

dyml <- yaml.load_file(args[1])
rds <- dyml$rds
meta <- dyml$meta
assay <- dyml$assay
ctl <- dyml$ctl
form <- dyml$form
ref_gene <- dyml$ref_gene
samples <- dyml$samples
outdir <- dyml$outdir

dir.create(outdir)


seurat_obj <- readRDS(rds)
metadata <- fread(meta, data.table=F)


cell_anno <- metadata[,c('cell_type2', 'orig.ident')]
cell_anno$cell <- metadata$V1

ref_group <- strsplit(ctl, ',')[[1]]
print(ref_group)


if (form == 'sample'){
    cell_anno$group <- cell_anno$orig.ident
}

if (form == 'celltype'){
    cell_anno$group <- cell_anno$cell_type2
}

cell_anno$group <- as.character(cell_anno$group)
cell_anno$group[cell_anno$cell_type2 %in% ref_group] <- 'Control'

# output cell anno data
cell_anno <- cell_anno[, c('cell', 'group')]
# extract by sample_group + all 'Countrol'
if (samples != 'all'){
    samples <- stringr::str_split(samples, ',')[[1]]
    cell_anno <- cell_anno %>% filter(cell_anno$group %in% c('Control', samples))
}
fwrite(cell_anno, file=paste0(outdir, '/cellInfo.txt'), sep='\t', row.names=F, col.names=F)




# exp count
data_genes <- row.names(seurat_obj@assays[[assay]]$counts)


# gene anno
d_gene <- read.table(ref_gene, sep='\t', header=F)
d_gene <- d_gene %>% distinct(V1, .keep_all = TRUE)
rownames(d_gene) <- d_gene$V1
d_gene$V1 <- NULL
used_genes <- intersect(data_genes, rownames(d_gene))
d_gene2 <- d_gene[used_genes,]
fwrite(d_gene2, file=paste0(outdir, '/geneLoci.txt'), sep = "\t", quote=F, row.names=T, col.names=F)


# exp count for used_genes
assay_count <- seurat_obj@assays[[assay]]$counts[used_genes, cell_anno$cell]
count_mtx <- as.data.frame(as.matrix(assay_count))
rownames(count_mtx) <- used_genes
fwrite(count_mtx, file=paste0(outdir, '/expCount.txt'), sep = "\t", quote=F, row.names=T, col.names=T)





