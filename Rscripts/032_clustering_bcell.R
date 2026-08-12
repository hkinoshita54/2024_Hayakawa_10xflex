# Continued from 030_clustering_imm.R
# Without integration
analysis_step <- "032_clustering_bcell"

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
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Load files ----
pal <- readRDS("RDSfiles/color_palette.RDS")
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load("RDSfiles/immgroup_names.Rdata")
seu <- seu_all[, B_names]
seu$cellgroup <- "Imm."
seu$celltype1 <- "Bcell"

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 0.5)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

# find markers, excluding obvious non-b genes
non_b_genes <- c("Krt18","Krt19","Muc4","Muc6","Tff2",
                 "Acta2","Col6a3","Sox9","Il1b","Lyz2","S100a8","S100a9")

all_genes   <- rownames(seu)
b_gene_pool <- setdiff(all_genes, non_b_genes)

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  features = b_gene_pool
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top30
write_csv(top30, file = file.path(res_path, "top30_markers.csv"))

# Check markers by feature plots
features <- readLines("aux_data/gene_set/additional/32_bcell_markers.txt")
sapply(features, save_fp, seu, fp_path)

add_feat <- "Tnfrsf17"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# adjust resolution
# seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters")
ggsave("cluster.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 45, height = 45, units = "mm")

# Add celltype annotation and save the Seurat object
seu$celltype2 <- ""
seu$celltype2[seu$seurat_clusters %in% c(0,1)] <- "Conv.-B"
seu$celltype2[seu$seurat_clusters %in% c(2)] <- "GC-B"
seu$celltype2[seu$seurat_clusters %in% c(3)] <- "Cycling-B"
seu$celltype2[seu$seurat_clusters %in% c(4)] <- "Plasma"
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "Pre-B"
seu$celltype2 <- factor(seu$celltype2, levels = c("Pre-B", "Conv.-B", "GC-B", "Cycling-B", "Plasma"))
DimPlot(seu, group.by = "celltype2", cols = pal$B_cells, raster = TRUE, raster.dpi = c(600, 600), pt.size = 6) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, width = 45, height = 45, units = "mm")

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  features = b_gene_pool
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top30
write_csv(top30, file = file.path(res_path, "top30_markers_annotated.csv"))
features <- c(
  "Rag1","Il7r",        # Pre-B
  "Grk6","Trim30b",     # Conv.-B
  "H2-Aa","Cxcr5",      # GC-B
  "Mki67","Top2a",      # Cycling-B
  "Jchain","Xbp1"       # Plasma
)
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
ggsave("dotplot.pdf", path = plot_path, width = 60, height = 45, units = "mm")

saveRDS(seu, file = file.path("RDSfiles", "seu_032_bcell.RDS"))

# Cell numbers and proportion plot ----
seu <- readRDS("RDSfiles/seu_032_bcell.RDS")
fraction <- "bcell"

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
  scale_fill_manual(values = pal$B_cells) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 40, height = 45, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$B_cells) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 45, units = "mm")

