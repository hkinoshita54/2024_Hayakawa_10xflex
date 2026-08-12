# Continued from 010.1_clustering_2.R
# Without integration 
analysis_step <- "030_clustering_imm"

# load packages ----
library(tidyverse)
library(Seurat)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper function ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
pal <- readRDS("RDSfiles/color_palette.RDS")
load(file = "RDSfiles/cellcycle_genes_mouse.Rdata")
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))
seu <- seu_all[, imm_names]
seu$cellgroup <- "Imm."

# Clustering ----
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- RunPCA(seu, npcs = 50)
seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, dims = 1:50)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top30
write_csv(top30, file = file.path(res_path, "top30_markers.csv"))
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# ggsave("QC_vln_unfil.png", path = plot_path, width = 15, height = 3, units = "in", dpi = 150)
# seu <- subset(seu, subset = seurat_clusters %in% c(10), invert = TRUE)

# adjust resolution
# seu <- FindClusters(seu, resolution = 3, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) & NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)


# Check markers by feature plots
features <- readLines(file.path("aux_data", "gene_set", "annotation", "02_imm_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- "Pgc"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype1 annotation
seu$celltype1 <- "Mye."
seu$celltype1[seu$seurat_clusters %in% c(5,6,18,19)] <- "Bcell"
seu$celltype1[seu$seurat_clusters %in% c(1,4,13,14)] <- "Tcell"
seu$celltype1 <- factor(seu$celltype1, levels = c("Bcell", "Tcell", "Mye."))
DimPlot(seu, group.by = "celltype1", label = TRUE, repel = TRUE, cols = pal$immune_group) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype1.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

# save the cell barcodes of each celltype1
B_names <- colnames(seu)[seu$celltype1 == "Bcell"]
T_names <- colnames(seu)[seu$celltype1 == "Tcell"]
Mye_names <- colnames(seu)[seu$celltype1 == "Mye."]
save(B_names, T_names, Mye_names, file = file.path("RDSfiles", "immgroup_names.Rdata"))

# Dot plot
Idents(seu) <- "celltype1"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Cd19", "Pax5", "Cd3e", "Cd3g", "Csf2ra", "Tyrobp")
DotPlot(seu, group.by = "celltype1", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 4.5, height = 4.5, units = "in", dpi = 150)
