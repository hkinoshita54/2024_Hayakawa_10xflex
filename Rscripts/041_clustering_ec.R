# Continued from 040_clustering_str.R
# Without integration
analysis_step <- "041_clustering_ec"

# load packages ----
library(tidyverse)
library(Seurat)
source("Rscripts/theme_panel.R")

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper function ----
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

pal <- readRDS("RDSfiles/color_palette.RDS")

# Load files ----
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load("RDSfiles/str_celltype1_names.Rdata")
seu <- seu_all[, EC_names]
seu$cellgroup <- "Str."
seu$celltype1 <- "EC"

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))
seu_cp <- seu

# remove contamination, then re-cluster
seu <- subset(seu, subset = seurat_clusters %in% 8, invert = TRUE) # contamination of epithelial cells
seu <- recluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_2.csv"))
seu <- subset(seu, subset = seurat_clusters %in% 8, invert = TRUE) # contamination of pericytes

# Check markers by feature plots
features <- readLines("aux_data/gene_set/additional/endothelial_markers.txt")
sapply(features, save_fp, seu, fp_path)

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
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "Art-EC"
seu$celltype2[seu$seurat_clusters %in% c(4)] <- "Ven-EC"
seu$celltype2[seu$seurat_clusters %in% c(0)] <- "Cap-EC"
seu$celltype2[seu$seurat_clusters %in% c(1)] <- "Tum-EC1"
seu$celltype2[seu$seurat_clusters %in% c(2)] <- "Tum-EC2"
seu$celltype2[seu$seurat_clusters %in% c(7)] <- "Tum-EC3"
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "LEC"
seu$celltype2[seu$seurat_clusters %in% c(3)] <- "Tum-LEC"
seu$celltype2 <- factor(seu$celltype2, levels = c("Art-EC", "Ven-EC", "Cap-EC", 
                                                  "Tum-EC1", "Tum-EC2", "Tum-EC3", "LEC", "Tum-LEC"))
DimPlot(seu, group.by = "celltype2", cols = pal$endothelial, raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

saveRDS(seu, file = file.path("RDSfiles", "seu_041_ec.RDS"))

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

## cluster markers
features <- c("Sema3g","Gja5",
              "Vwf","Eln",
              "Cd300lg","Car4",
              "Kit","Serpine1",
              "Selp","Glycam1",
              "Top2a","Birc5",
              "Prox1","Lyve1",
              "Flt4","Ackr2")
features <- c("Sema3g",
              "Vwf",
              "Cd300lg",
              "Serpine1",
              "Glycam1",
              "Top2a",
              "Lyve1")
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  ) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot2.pdf", path = plot_path, width = 45, height = 36, units = "mm")

## tan accepting
adhesion <- c("Icam1", "Vcam1", "Selp", "Sele")
# chemokines <- c("Csf3", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5")
secreted <- c("Pdgfa","Pdgfb", "Sema3f")
DotPlot(seu, group.by = "celltype2", features = unique(c("Eng", "Mcam", "Cdh5", adhesion, secreted)), dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  ) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot_sel.pdf", path = plot_path, width = 60, height = 35, units = "mm")


# EDA ----
seu <- readRDS(file = file.path("RDSfiles", "seu_041_ec.RDS"))

add_feat <- c("Col3a1")
sapply(add_feat, save_fp, seu, fp_path)

VlnPlot(seu, group.by = "celltype2", features = add_feat, cols = pal$endothelial, pt.size = 0) &
  theme_panel() & NoLegend() &
  labs(x = NULL, y = "Expression") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) 
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 35, height = 35, units = "mm")

## check the distribution of key CCC genes in seurat object
features <- c("Lrg1", "Eng", "Grem1", "Wnt5a", "Tgfb1", 
              "Spp1", "Mcam", "Cdh5",
              "Pdgfra", "Itgb4", "Itga6",
              "Plau", "Mdk", "Angptl4", "Lgals9")
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
ggsave("dotplot_ccc_genes.pdf", path = plot_path, width = 75, height = 50, units = "mm")




# Cell numbers and proportion plot ----
library(speckle)
fraction <- "ec"

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
  scale_fill_manual(values = pal$endothelial) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 42, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$endothelial) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 28, height = 35, units = "mm")
