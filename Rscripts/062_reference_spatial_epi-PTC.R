# doing this AFTER 025
analysis_step <- "062_reference_spatial_epi-PTC"

# load packages ----
library(tidyverse)
library(Seurat)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

# load data ----
ptc_cols <- c(PTC1 = "#CC79A7", PTC2 = "#7F7F7F", PTC3 = "#E69F00", PTC4 = "#009E73")
seu <- readRDS("RDSfiles/seu_025_ptc.RDS")
DimPlot(seu, group.by = "celltype3.1", cols = ptc_cols) & NoAxes()

# select genes for rctd ----
## 2000 hvg + top 50 marker genes for each cluster - mito, ribo
seu <- FindVariableFeatures(seu, nfeatures = 2000)
hvg <- VariableFeatures(seu)

Idents(seu) <- "celltype3.1"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers_use <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  pull(gene) %>% 
  unique()

genes_use <- unique(c(hvg, markers_use))
genes_use <- genes_use[!grepl("^mt-|^Rpl|^Rps", genes_use, ignore.case=FALSE)]
# load("RDSfiles/cellcycle_genes_mouse.Rdata")
# genes_use <- setdiff(genes_use, unique(c(s.genes, g2m.genes)))

## save for later use as a reference
metadata <- seu@meta.data %>% 
  select("orig.ident", "cellgroup", "celltype1", "celltype2", "celltype3.1")
cts <- LayerData(seu, assay = 'RNA', layer = 'counts')
seu_ref <- CreateSeuratObject(cts, meta.data = metadata)
saveRDS(seu_ref, "RDSfiles/seu_062_ref_epi-PTC.RDS")

seu_rctd <- seu_ref[genes_use,]
cts_rctd <- LayerData(seu_rctd, assay = 'RNA', layer = 'counts')
seu_rctd <- CreateSeuratObject(cts_rctd, meta.data = metadata)
saveRDS(seu_rctd, "RDSfiles/seu_062_ref_epi-PTC_rctd.RDS")
