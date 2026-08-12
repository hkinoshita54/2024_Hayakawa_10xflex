# Modify output directory and seurat object to read in

# load packages ----
library(tidyverse)
library(Seurat)

# output directory
description <- "053_combine_epi"
out_dir <- file.path("out_data", description)
fs::dir_create(out_dir)

# load data
seu <- readRDS("RDSfiles/seu_053_epi_combined.RDS")
DimPlot(seu, group.by = "celltype2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()

# convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu$barcode <- colnames(seu)
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
seu <- FindVariableFeatures(seu, nfeatures = 8000)
hvgs <- VariableFeatures(seu)
tfs <- readLines("pyscenic_docker/aux_data/allTFs_mm.txt")
tfs <- tfs[tfs %in% rownames(seu)]
genes_for_scenic <- union(hvgs, tfs)
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
counts_matrix <- counts_matrix[genes_for_scenic,]
# counts_matrix <- counts_matrix[hvgs,]
writeMM(counts_matrix, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)
