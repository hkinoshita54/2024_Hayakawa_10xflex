# Continued from 030_clustering_imm.R
# Without integration
analysis_step <- "031_clustering_tcell"

# load packages ----
library(tidyverse)
library(Seurat)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper function ----
source("Rscripts/theme_panel.R")
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

# Load files ----
pal <- readRDS("RDSfiles/color_palette.RDS")
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load("RDSfiles/immgroup_names.Rdata")
seu <- seu_all[, T_names]
seu$cellgroup <- "Imm."
seu$celltype1 <- "Tcell"

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
markers <- FindAllMarkers(seu, only.pos = TRUE) # c1 seems to be epithelial contamination > leave as it is
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top30
write_csv(top30, file = file.path(res_path, "top30_markers.csv"))

# Check markers by feature plots
features <- readLines("aux_data/gene_set/additional/31_tnkcell_markers.txt")
sapply(features, save_fp, seu, fp_path)

add_feat <- "Rora"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# adjust resolution
# seu <- FindClusters(seu, resolution = 1.5, verbose = FALSE)
DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters")
ggsave("cluster.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 45, height = 45, units = "mm")

# Add celltype annotation and save the Seurat object
seu$celltype2 <- ""
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "Cycling-T"
seu$celltype2[seu$seurat_clusters %in% c(3)] <- "Stem-T"
seu$celltype2[seu$seurat_clusters %in% c(2)] <- "CD4-T"
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "CD8-T"
seu$celltype2[seu$seurat_clusters %in% c(8)] <- "Treg"
seu$celltype2[seu$seurat_clusters %in% c(9)] <- "NK"
seu$celltype2[seu$seurat_clusters %in% c(1)] <- "ILC2"
seu$celltype2[seu$seurat_clusters %in% c(0,4,7)] <- "gdT"
seu$celltype2 <- factor(seu$celltype2, levels = c("Cycling-T", "Stem-T", "CD4-T", "CD8-T", "Treg", "NK", "ILC2", "gdT"))
DimPlot(seu, group.by = "celltype2", cols = pal$T_cells, raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, width = 35, height = 30, units = "mm")

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Top2a", "Mki67", "Tcf7", "Lef1", "Cd4", "Ccr7", "Cd8a", "Cxcr3",
              "Foxp3", "Ctla4", "Ncr1", "Klrb1c", "Gata3", "Il1rl1", "Trdc", "Sox13")
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -2, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 70, height = 45, units = "mm")

saveRDS(seu, file = file.path("RDSfiles", "seu_031_tcell.RDS"))
seu <- readRDS(file = file.path("RDSfiles", "seu_031_tcell.RDS"))

# Cell numbers and proportion plot ----
fraction <- "tcell"

nclust <- length(levels(seu$celltype2))
props <- getTransformedProps(clusters = seu$celltype2, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$T_cells) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 40, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$T_cells) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 35, units = "mm")


# EDA ----
add_feat <- c("Il17a", "Il17f", "Rorc", "Ccr6", "Ifng", "Tbx21", "Nkg7", "Gzmb")
sapply(add_feat, save_fp, seu, fp_path)
