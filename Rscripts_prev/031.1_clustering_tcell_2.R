# Continued from 030.1_clustering_imm_2.R
# Without integration
analysis_step <- "031.1_clustering_tcell_2"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper function ----
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

# Clustering ----
seu <- readRDS(file.path("RDSfiles", "seu_030.1_imm.RDS"))
seu <- subset(seu, subset = immgroup == "T&NK")
seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
markers <- FindAllMarkers(seu, only.pos = TRUE) # c1 seems to be epithelial contamination > leave as it is

# adjust resolution
seu <- FindClusters(seu, resolution = 1.5, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)
markers <- FindAllMarkers(seu, only.pos = TRUE)

DimPlot(seu, label = TRUE, repel = TRUE, group.by = "celltype", cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

features <- c("Sell", "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Ctla4", "Trdc", "Rorc", 
              "Ncr1",  "Nkg7", "Top2a", "Mki67", "Gata3", "Rora")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot_celltype.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)


# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "31_tnkcell_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- "Rora"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype annotation and save the Seurat object
seu$celltype2 <- ""
seu$celltype2[seu$seurat_clusters %in% c(2)] <- "CD4_Tn"
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "CD4-T"
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "CD8-T"
seu$celltype2[seu$seurat_clusters %in% c(7)] <- "Treg"
seu$celltype2[seu$seurat_clusters %in% c(0,3,8,9)] <- "gdT"
seu$celltype2[seu$seurat_clusters %in% c(10)] <- "NK"
seu$celltype2[seu$seurat_clusters %in% c(1)] <- "ILC2"
seu$celltype2[seu$seurat_clusters %in% c(4)] <- "Prolif.T"
seu$celltype2 <- factor(seu$celltype2, levels = c("CD4_Tn", "CD4-T", "CD8-T", "Treg", "gdT", "NK", "ILC2", "Prolif.T"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", label = TRUE, repel = TRUE) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype2.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Ccr7", "Lef1", "Sell", "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Ctla4", "Trdc", "Rorc", 
              "Ncr1",  "Nkg7", "Gata3", "Rora", "Top2a", "Mki67")
DotPlot(seu, group.by = "celltype2", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_031.1_tcell.RDS"))
