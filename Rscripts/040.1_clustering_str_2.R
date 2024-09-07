# Continued from 010.1_clustering_2.R
analysis_step <- "040.1_clustering_str_2"

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
seu <- seu_all[, str_names]
seu$cellgroup <- "Str."

seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
markers <- FindAllMarkers(seu, only.pos = TRUE) # Check markers

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "03_str_markers.txt"))
sapply(features, save_fp, seu, fp_path)

features <- readLines(file.path("aux_data", "gene_set", "annotation", "endothelial_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Cxcl12")
sapply(add_feat, save_fp, seu, fp_path)

# annotation ----
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(0,2,3,5,6,8,12,16,17,18,24,26,29)] <- "Fib."
seu$celltype[seu$seurat_clusters %in% c(4,7,13,15,20,27)] <- "Myo."
seu$celltype[seu$seurat_clusters %in% c(1,9,11,19,21)] <- "Endo."
seu$celltype[seu$seurat_clusters %in% c(10)] <- "Peri."
seu$celltype[seu$seurat_clusters %in% c(14)] <- "Glia."
seu$celltype[seu$seurat_clusters %in% c(22)] <- "Adipo."
seu$celltype[seu$seurat_clusters %in% c(23)] <- "Meso."
seu$celltype[seu$seurat_clusters %in% c(25)] <- "ICC"
seu$celltype[seu$seurat_clusters %in% c(28)] <- "Neuro."
seu$celltype <- factor(seu$celltype, levels = c("Fib.", "Myo.", "Endo.", "Peri.", "Glia.", "Adipo.", "Meso.", "ICC", "Neuro."))
DimPlot(seu, group.by = "celltype", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1000) %>%
  ungroup() -> markers
openxlsx2::write_xlsx(markers, file.path(res_path, "markers.xlsx"))
features <- c("Pdgfra", "Myh11", "Pecam1", "Rgs5", "Sox10", "Plin1", "Msln", "Kit", "Snap25")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_040.1_str.RDS"))
