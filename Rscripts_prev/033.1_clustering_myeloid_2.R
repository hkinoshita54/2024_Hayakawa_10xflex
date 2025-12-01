# Continued from 030.1_clustering_imm_2.R
# Without integration
analysis_step <- "033.1_clustering_myeloid_2"

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
seu <- readRDS(file.path("RDSfiles", "seu_030.2_imm.RDS"))
seu <- subset(seu, subset = immgroup == "Myeloids")
seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
markers <- FindAllMarkers(seu, only.pos = TRUE) # There seems to be epithelial contamination > leave as it is

# adjust resolution
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
# markers <- FindAllMarkers(seu, only.pos = TRUE)

DimPlot(seu, label = TRUE, repel = TRUE, group.by = "celltype", cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

features <- c("C1qc", "Apoe", "Thbs1", "Inhba", "Vcan", "F10",  
              "Clec9a",  "Batf3", "Siglech", "Csf3r", "S100a9","Tpsb2", "Ms4a2")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot_celltype.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "33_myeloid_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- "Col4a2"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype annotation and save the Seurat object
seu$celltype2 <- ""
seu$celltype2[seu$seurat_clusters %in% c(1,7,24,25,26)] <- "Mac-1"
seu$celltype2[seu$seurat_clusters %in% c(20)] <- "Mac-2" #resident
seu$celltype2[seu$seurat_clusters %in% c(3,9,17,22,23)] <- "TAM-1" #Ctsd
seu$celltype2[seu$seurat_clusters %in% c(11)] <- "TAM-2" #Cd163
seu$celltype2[seu$seurat_clusters %in% c(2)] <- "TAM-3" #Inhba
seu$celltype2[seu$seurat_clusters %in% c(4,10,15)] <- "Mono-1" #Vcan
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "Mono-2" #Other
seu$celltype2[seu$seurat_clusters %in% c(0,6,8,13,18)] <- "TAN"
seu$celltype2[seu$seurat_clusters %in% c(12,14,19)] <- "cDC"
seu$celltype2[seu$seurat_clusters %in% c(21)] <- "pDC"
seu$celltype2[seu$seurat_clusters %in% c(16)] <- "Mast"
seu$celltype2 <- factor(seu$celltype2, levels = c("Mac-1", "Mac-2", "TAM-1", "TAM-2", "TAM-3", 
                                                  "Mono-1", "Mono-2", "TAN", "cDC", "pDC", "Mast"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", label = TRUE, repel = TRUE) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype2.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("C1qc", "Cd163", "Apoe", "Ctsd", "Thbs1", "Inhba", "Vcan", "F10", "Csf3r", "S100a9", 
              "Clec9a",  "Batf3", "Siglech", "Tpsb2", "Ms4a2")
DotPlot(seu, group.by = "celltype2", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_033.1_myeloid.RDS"))
