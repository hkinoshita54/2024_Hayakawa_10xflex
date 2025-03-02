analysis_step <- "910_reference_for_cosmx"

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

# Clustering ----
seu <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))

seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
# ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150) 
DimPlot(seu, group.by = "orig.ident") + NoAxes()
# ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
# markers <- FindAllMarkers(seu, only.pos = TRUE)

# fp_path <- file.path(plot_path, "feature_plot")
# fs::dir_create(c(fp_path))
# files <- list.files(path = file.path("aux_data", "gene_set", "annotation"), pattern = ".txt", full.names = TRUE)
# features <- lapply(files, FUN = readLines) %>% unlist()
# sapply(features, save_fp, seu, fp_path)

# Add cellgroup annotation ----
epi_names <- colnames(seu)[seu$seurat_clusters %in% c(0:4,17,19,20,22,25,27,29,30,31,34)]
imm_names <- colnames(seu)[seu$seurat_clusters %in% c(7,9,10,13,16,18,21,33)]
str_names <- colnames(seu)[seu$seurat_clusters %in% c(5,6,8,11,12,14,15,23,24,26,28,32,35)]
# save(epi_names, str_names, imm_names, file = file.path("RDSfiles", "cellgroup_names_2.Rdata"))

seu$cellgroup <- "Epi."
seu$cellgroup[imm_names] <- "Imm."
seu$cellgroup[str_names] <- "Str."
seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Imm.", "Str."))
DimPlot(seu, group.by = "cellgroup") & NoAxes()
# ggsave("cellgroup.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
# Idents(seu) <- "cellgroup"
# markers <- FindAllMarkers(seu, only.pos = TRUE)
# features <- c("Epcam", "Cldn18", "Ptprc", "Arhgap45", "Dcn", "Sparc")
# DotPlot(seu, group.by = "cellgroup", features = features) + RotatedAxis()
# ggsave("dotplot.png", path = plot_path, width = 4, height = 3.7, units = "in", dpi = 150)

seu <- subset(seu, subset = orig.ident %in% c("WT", "PTC"))
DimPlot(seu, group.by = "cellgroup") & NoAxes()
DimPlot(seu, group.by = "orig.ident") + NoAxes()
saveRDS(seu, file = file.path("RDSfiles", "seu_910_reference_for_cosmx.RDS"))
