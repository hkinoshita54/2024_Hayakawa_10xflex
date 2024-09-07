# Continued from 010.1_clustering_2.R
# Without integration
analysis_step <- "030.1_clustering_imm_2"

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
seu <- seu_all[, imm_names]
seu$cellgroup <- "Imm."

seu <- cluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE) # no apparent contamination
# seu <- subset(seu, subset = seurat_clusters %in% 20:22, invert = TRUE)
# seu <- recluster(seu, npcs = 50, res = 2)

seu <- FindClusters(seu, resolution = 3, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))

ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "02_imm_markers.txt"))
sapply(features, save_fp, seu, fp_path)

features <- readLines(file.path("aux_data", "gene_set", "annotation", "additional_imm_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- "Cxcl12"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype annotation and save the Seurat object
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(2,5,10,35)] <- "Bcell"
seu$celltype[seu$seurat_clusters %in% c(31,40)] <- "Plasma"
seu$celltype[seu$seurat_clusters %in% c(18,21)] <- "CD4-T"
seu$celltype[seu$seurat_clusters %in% c(22)] <- "CD8-T"
seu$celltype[seu$seurat_clusters %in% c(23)] <- "Treg"
seu$celltype[seu$seurat_clusters %in% c(8,17,28)] <- "gdT"
seu$celltype[seu$seurat_clusters %in% c(36)] <- "NK"
seu$celltype[seu$seurat_clusters %in% c(19)] <- "Prolif.T"
seu$celltype[seu$seurat_clusters %in% c(20)] <- "ILC2"
seu$celltype[seu$seurat_clusters %in% c(4,11,32,37,38,39)] <- "Macro."
seu$celltype[seu$seurat_clusters %in% c(3,14,16,29)] <- "TAM"
seu$celltype[seu$seurat_clusters %in% c(1,6,7,13,26)] <- "Mono."
seu$celltype[seu$seurat_clusters %in% c(15,25,30)] <- "cDC"
seu$celltype[seu$seurat_clusters %in% c(34)] <- "pDC"
seu$celltype[seu$seurat_clusters %in% c(0,9,12,24,33)] <- "TAN"
seu$celltype[seu$seurat_clusters %in% c(27)] <- "Mast"
seu$celltype <- factor(seu$celltype, levels = c("Bcell", "Plasma", "CD4-T", "CD8-T", "Treg", "gdT", "NK", "Prolif.T",
                                                "ILC2", "Macro.", "TAM", "Mono.", "cDC", "pDC", "TAN", "Mast"))
DimPlot(seu, group.by = "celltype", cols = "polychrome", label = TRUE, repel = TRUE) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Cd79a", "Jchain", "Cd4", "Cd8a", "Foxp3", "Trdc", "Nkg7", "Top2a",
              "Gata3", "C1qb", "Il1b", "Batf3", "Siglech", "S100a8", "Kit")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_030.1_imm.RDS"))
