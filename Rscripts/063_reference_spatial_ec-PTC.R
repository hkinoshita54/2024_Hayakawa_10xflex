# doing this AFTER 041
analysis_step <- "063_reference_spatial_ec-PTC"

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
pal <- readRDS("RDSfiles/color_palette.RDS")
seu_all <- readRDS("RDSfiles/seu_041_ec.RDS")
table(seu_all$orig.ident, seu_all$celltype2)

# subset to PTC, with the exception for Squam., Bcell, Glial and Neural ----
## drop Stem, Pit, EEC, Tuft, Adipo., Meso. which we don't analyze later
seu <- subset(seu_all, subset = orig.ident == "PTC")
seu <- subset(seu, subset = celltype2 == "Art-EC", invert = TRUE)
seu$celltype2 <- droplevels(seu$celltype2)
table(seu$celltype2)

# plots for record ----
seu <- cluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) + NoAxes()
ggsave("sample.png", path = plot_path, width = 3.3, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype2", cols = "polychrome") & NoAxes()
ggsave("celltype2.png", path = plot_path, width = 3.4, height = 3, units = "in", dpi = 150)

# select genes for rctd ----
## 2000 hvg + top 50 marker genes for each cluster - mito, ribo, cell cycle
hvg <- VariableFeatures(seu)

Idents(seu) <- "celltype2"
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
  select("orig.ident", "cellgroup", "celltype1", "celltype2")
cts <- LayerData(seu, assay = 'RNA', layer = 'counts')
seu_ref <- CreateSeuratObject(cts, meta.data = metadata)
saveRDS(seu_ref, "RDSfiles/seu_063_ref_ec-PTC.RDS")

seu_rctd <- seu_ref[genes_use,]
cts_rctd <- LayerData(seu_rctd, assay = 'RNA', layer = 'counts')
seu_rctd <- CreateSeuratObject(cts_rctd, meta.data = metadata)
saveRDS(seu_rctd, "RDSfiles/seu_063_ref_ec-PTC_rctd.RDS")
