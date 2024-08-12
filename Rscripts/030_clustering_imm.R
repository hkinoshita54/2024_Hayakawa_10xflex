# Continued from 010_QC.R
# Without integration
analysis_step <- "030_clustering_imm"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))

# Clustering w/o integration ----
plot_path <- file.path("plot", analysis_step, "no_int")    # set a new plot directory
fs::dir_create(plot_path)
seu <- seu_all[, imm_names]
seu$cellgroup <- "Imm."
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)    # adjust resolution when necessary
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "annotation", "02_imm_markers.txt"))

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path,
           width = 4, height = 4, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
# markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Tcl1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

seu_wo_int <- seu

# Clustering w/ harmony integration ----
plot_path <- file.path("plot", analysis_step, "harmony")    # set a new plot directory
fs::dir_create(plot_path)
seu <- seu_all[, imm_names]
seu$cellgroup <- "Str."
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident) 
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
# seu <- JoinLayers(seu)
# markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Bmp5"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)

seu_harm <- seu

# Clustering w/ Harmony integration only for batch (WT, PT, PC as 2nd) ----
plot_path <- file.path("plot", analysis_step, "harmony_batch")    # set a new plot directory
fs::dir_create(plot_path)
seu <- seu_all[, imm_names]
seu$cellgroup <- "Str."
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$batch) 
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
# seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Jchain"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)

seu_harm_batch <- seu

# Add celltype annotation and save the Seurat object ----
# use clustering w/o integration
plot_path <- file.path("plot", analysis_step, "no_int")    # set a new plot directory
seu <- seu_wo_int
seu <- subset(seu, idents = c(7), invert = TRUE)    # These clusters seem contamination
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(0,6)] <- "Bcell"
seu$celltype[seu$seurat_clusters %in% c(27)] <- "Plasma"
seu$celltype[seu$seurat_clusters %in% c(5)] <- "CD4-T"
seu$celltype[seu$seurat_clusters %in% c(18)] <- "CD8-T"
seu$celltype[seu$seurat_clusters %in% c(21)] <- "Treg"
seu$celltype[seu$seurat_clusters %in% c(1,22)] <- "gdT"
seu$celltype[seu$seurat_clusters %in% c(19)] <- "ILC2"
seu$celltype[seu$seurat_clusters %in% c(15)] <- "Prolif.T&B"
seu$celltype[seu$seurat_clusters %in% c(8,9,13,14,16,24,25)] <- "Macro."
seu$celltype[seu$seurat_clusters %in% c(3,10,11)] <- "Infl.Mono."
seu$celltype[seu$seurat_clusters %in% c(12,20)] <- "cDC"
seu$celltype[seu$seurat_clusters %in% c(28)] <- "pDC"
seu$celltype[seu$seurat_clusters %in% c(2,4,17,26)] <- "Neutro."
seu$celltype[seu$seurat_clusters %in% c(23)] <- "Mast"
seu$celltype <- factor(seu$celltype, levels = c("Bcell", "Plasma", "CD4-T", "CD8-T", "Treg", "gdT",  
                                                "ILC2", "Prolif.T&B", "Macro.", "Infl.Mono.", "cDC", "pDC", "Neutro.", "Mast"))
DimPlot(seu, group.by = "celltype", cols = "alphabet", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Cd79a", "Jchain", "Cd4", "Cd8a", "Foxp3", "Trdc", 
              "Gata3", "Top2a", "C1qb", "Vcan", "Il1b", "Clec9a", "Siglech", "S100a8", "Kit")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6, height = 4, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_030_imm.RDS"))

