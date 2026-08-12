# Continued from 021_clustering_epi-PC&PTC.R
# 2026-01-10
analysis_step <- "025_clustering_epi-PTC"

# load packages ----
library(tidyverse)
library(Seurat)
library(speckle)
library(openxlsx2)
source("Rscripts/theme_panel.R")

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 5) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = fp_path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# set color palettes ----
pal <- readRDS("RDSfiles/color_palette.RDS")
ptc_cols <- c(PTC1 = "#CC79A7", PTC2 = "#7F7F7F", PTC3 = "#E69F00", PTC4 = "#009E73")

# Clustering ----

## start from 021, whose contamination are already removed
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
seu <- subset(seu, subset = celltype2 == "Epi-PTC" & orig.ident == "PTC")
seu$orig.ident <- droplevels(seu$orig.ident)
seu <- recluster(seu, npcs = 50, res = 0.3)

DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Cluster")
ggsave("cluster.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "celltype3", cols = "polychrome", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype3")
ggsave("celltype3.pdf", path = plot_path, width = 45, height = 45, units = "mm")
table(seu$celltype3)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.3.csv"))

## Add celltyp3.1 annotation and save the Seurat object ----
seu$celltype3.1 <- ""
seu$celltype3.1[seu$seurat_clusters %in% c(2)] <- "PTC1" # secretory & proliferative
seu$celltype3.1[seu$seurat_clusters %in% c(0)] <- "PTC2" # NOS
seu$celltype3.1[seu$seurat_clusters %in% c(1)] <- "PTC3" # invasive
seu$celltype3.1[seu$seurat_clusters %in% c(3)] <- "PTC4" # inflammatory
seu$celltype3.1 <- factor(seu$celltype3.1, levels = c("PTC1", "PTC2", "PTC3", "PTC4"))
DimPlot(seu, group.by = "celltype3.1", cols = ptc_cols, raster = TRUE, pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "celltype3.1")
ggsave("celltype3.1.pdf", path = plot_path, width = 30, height = 30, units = "mm")

Idents(seu) <- "celltype3.1"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))
write_csv(markers, file = file.path(res_path, "markers_annotated.csv"))
openxlsx2::write_xlsx(markers, file = file.path(res_path, "markers_annotated.xlsx"))

saveRDS(seu, file = "RDSfiles/seu_025_ptc.RDS")

# EDA ----
seu <- readRDS("RDSfiles/seu_025_ptc.RDS")

# add_feat <- c("Lrg1", "Cd38", "Cldn4", "Msln", "Cxcl5", "Pigr")
add_feat <- c("Spp1")
sapply(add_feat, save_fp, seu, fp_path)

VlnPlot(seu, group.by = "celltype3.1", features = add_feat, cols = ptc_cols, pt.size = 0) &
  theme_panel() & NoLegend() &
  labs(x = NULL, y = "Expression") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 25, height = 30, units = "mm")

add_feat <- "Nt5e"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 5) &
  theme_panel() & NoAxes() & NoLegend()
ggsave(paste0(add_feat, ".pdf"), path = fp_path, width = 25, height = 30, units = "mm")

features <- c("Cldn4", "Nt5e", "Mki67", "Top2a")
DotPlot(seu, group.by = "celltype3.1", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -2, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot2.pdf", path = plot_path, width = 45, height = 40, units = "mm")

## check the distribution of key CCC genes in seurat object
features <- c("Lrg1", "Eng", "Grem1", "Wnt5a", "Tgfb1", 
              "Spp1", "Mcam", "Cdh5",
              "Pdgfra", "Itgb4", "Itga6",
              "Plau", "Mdk", "Angptl4", "Lgals9")
DotPlot(seu, group.by = "celltype3.1", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -2, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot_ccc_genes.pdf", path = plot_path, width = 70, height = 45, units = "mm")



# Cell numbers and proportion plot ----
# seu <- readRDS("RDSfiles/seu_025_ptc.RDS")
fraction <- "ptc"

nclust <- length(levels(seu$celltype3.1))
props <- getTransformedProps(clusters = seu$celltype3.1, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = ptc_cols) +
  theme_classic()
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 1.85, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = ptc_cols) +
  theme_classic()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 1.9, height = 3.2, units = "in", dpi = 150)
