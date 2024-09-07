# Continued from 010.1_clustering_2.R
# Without integration
analysis_step <- "020.1_clustering_epi_2"

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
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names_2.Rdata"))
seu <- seu_all[, epi_names]
seu$cellgroup <- "Epi."

seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE) # c21, c22 - unknown, c20 - enterocyte
seu <- subset(seu, subset = seurat_clusters %in% 20:22, invert = TRUE)

seu <- recluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "01_epi_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Lrg1", "Cd38", "Cldn4", "Msln", "Cxcl5", "Pigr")
sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE)

seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res2.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

add_feat <- "Chga"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype annotation and save the Seurat object
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(13)] <- "Prog."
seu$celltype[seu$seurat_clusters %in% c(4,7,10)] <- "Pit"
seu$celltype[seu$seurat_clusters %in% c(0,2)] <- "Neck"
seu$celltype[seu$seurat_clusters %in% c(6,12,14,16)] <- "Chief"
seu$celltype[seu$seurat_clusters %in% c(8,18,19)] <- "Pariet."
seu$celltype[seu$seurat_clusters %in% c(27,28,29)] <- "EEC"
seu$celltype[seu$seurat_clusters %in% c(30)] <- "Tuft"
seu$celltype[seu$seurat_clusters %in% c(20,24)] <- "Squam."
seu$celltype[seu$seurat_clusters %in% c(22,25)] <- "Epi-PT"
seu$celltype[seu$seurat_clusters %in% c(1,3,5,9,11,15,17,21,23,26)] <- "Epi-PC&PTC"
seu$celltype <- factor(seu$celltype, levels = c("Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                "Squam.", "Epi-PT", "Epi-PC&PTC"))
DimPlot(seu, group.by = "celltype", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Mki67", "Muc5ac", "Muc6", "Cblif", "Atp4b", "Chga", "Dclk1", "Krt5", "Ces2a", "Pigr")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_020.1_epi.RDS"))
