analysis_step <- "800_PGE2_related_genes"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper functions ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}
recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


# Epithelial ----
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() + NoLegend()
DimPlot(seu, group.by = "orig.ident") + NoAxes()

fp_path <- file.path(plot_path, "feature_plot/epi")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/additional/PGE_related.txt")
sapply(features, save_fp, seu, fp_path)

dp_path <- file.path(plot_path, "dot_plot/epi")
fs::dir_create(c(dp_path))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = dp_path, width = 7, height = 4, units = "in", dpi = 150)

vp_path <- file.path(plot_path, "vln_plot/epi")
fs::dir_create(c(vp_path))
VlnPlot(seu, features = features, stack = T, flip = T) & NoLegend()
ggsave("vlnplot.png", path = vp_path, width = 4, height = 6, units = "in", dpi = 150)


# Stromal ----
seu <- readRDS(file.path("RDSfiles", "seu_040.1_str.RDS"))
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() + NoLegend()
DimPlot(seu, group.by = "orig.ident") + NoAxes()

fp_path <- file.path(plot_path, "feature_plot/str")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/additional/PGE_related.txt")
sapply(features, save_fp, seu, fp_path)

dp_path <- file.path(plot_path, "dot_plot/str")
fs::dir_create(c(dp_path))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = dp_path, width = 7, height = 4, units = "in", dpi = 150)

vp_path <- file.path(plot_path, "vln_plot/str")
fs::dir_create(c(vp_path))
VlnPlot(seu, features = features, stack = T, flip = T) & NoLegend()
ggsave("vlnplot.png", path = vp_path, width = 4, height = 6, units = "in", dpi = 150)

# Immune ----
seu <- readRDS(file.path("RDSfiles", "seu_030.1_imm.RDS"))
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() + NoLegend()
DimPlot(seu, group.by = "orig.ident") + NoAxes()

fp_path <- file.path(plot_path, "feature_plot/imm")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/additional/PGE_related.txt")
sapply(features, save_fp, seu, fp_path)

dp_path <- file.path(plot_path, "dot_plot/imm")
fs::dir_create(c(dp_path))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = dp_path, width = 7, height = 6, units = "in", dpi = 150)

vp_path <- file.path(plot_path, "vln_plot/imm")
fs::dir_create(c(vp_path))
VlnPlot(seu, features = features[-8], stack = T, flip = T) & NoLegend()
ggsave("vlnplot.png", path = vp_path, width = 6, height = 6, units = "in", dpi = 150)
