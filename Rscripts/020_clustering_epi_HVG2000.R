# Continued from 010.1_clustering.R
# Without integration
analysis_step <- "020_clustering_epi"

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

recluster_sct = function(seu_obj, npcs, res, n.neighbors, min.dist){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- SCTransform(seu, vst.flavor = "v2", verbose = F)
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs, n.neighbors = n.neighbors, min.dist = min.dist)
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
load(file.path("RDSfiles", "cellgroup_names.Rdata"))
seu <- seu_all[, epi_names]
seu$cellgroup <- "Epi."

seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE) # c24, c25 - unknown, c23 - enterocyte
seu <- subset(seu, subset = seurat_clusters %in% 23:25, invert = TRUE)

seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE)
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## see how they are clustered with SCTransform
# seu <- recluster_sct(seu, npcs = 20, res = 2, n.neighbors = 50, min.dist = 0.4)
# DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
# ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
# DimPlot(seu, group.by = "orig.ident") + NoAxes()
# ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "01_epi_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Lrg1", "Cd38", "Cldn4", "Msln", "Cxcl5", "Pigr")
sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE)

seu <- recluster(seu, npcs = 50, res = 3)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res2.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

add_feat <- "Il6"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype1 annotation
seu$celltype1 <- ""
seu$celltype1[seu$seurat_clusters %in% c(18)] <- "Stem"
seu$celltype1[seu$seurat_clusters %in% c(6,31)] <- "Prog."
seu$celltype1[seu$seurat_clusters %in% c(10,14,15)] <- "Pit"
seu$celltype1[seu$seurat_clusters %in% c(0,1)] <- "Neck"
seu$celltype1[seu$seurat_clusters %in% c(5,11,17,19)] <- "Chief"
seu$celltype1[seu$seurat_clusters %in% c(9,13,22)] <- "Pariet."
seu$celltype1[seu$seurat_clusters %in% c(26,28,30)] <- "EEC"
seu$celltype1[seu$seurat_clusters %in% c(32)] <- "Tuft"
seu$celltype1[seu$seurat_clusters %in% c(21,24)] <- "Squam."
seu$celltype1[seu$seurat_clusters %in% c(23,25)] <- "Pit-PT"
seu$celltype1[seu$seurat_clusters %in% c(2,3,4,7,8,12,16,20,27,29)] <- "Epi-PC&PTC"
seu$celltype1 <- factor(seu$celltype1, levels = c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                "Squam.", "Pit-PT", "Epi-PC&PTC"))
DimPlot(seu, group.by = "celltype1", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()

# Add celltype2 annotation
tum_names <- colnames(seu)[seu$celltype1 %in% c("Epi-PC&PTC")]
seu$celltype2 <- as.character(seu$celltype1)
seu$celltype2[tum_names] <- paste0("Epi-", seu$orig.ident[tum_names])
seu <- subset(seu, subset = celltype2 %in% c("Epi-WT", "Epi-PT"), invert = TRUE)
seu$celltype2 <- factor(seu$celltype2, levels = c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                  "Squam.", "Pit-PT", "Epi-PC", "Epi-PTC"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype1", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype1.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype1"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("markers_epi_celltype1.xlsx")))
features <- c("Mki67", "Muc5ac", "Muc6", "Cblif", "Atp4b", "Chga", "Dclk1", "Krt5", "Ces2a", "Pigr")
DotPlot(seu, group.by = "celltype1", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

seu$RNA_snn_res.1 <- NULL
saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))
