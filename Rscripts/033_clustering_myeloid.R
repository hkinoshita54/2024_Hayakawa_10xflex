# Continued from 030_clustering_imm.R
# Without integration
analysis_step <- "033_clustering_myeloid"

# load packages ----
library(tidyverse)
library(Seurat)
library(UCell)

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
seu <- seu_all[, Mye_names]
seu$cellgroup <- "Imm."
seu$celltype1 <- "Mye."

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))

# find markers with all genes
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top30
write_csv(top30, file = file.path(res_path, "top30_markers.csv"))

# find markers, excluding obvious non-mye genes
epi.genes <- c("Atp4a","Atp4b","Pgc","Clps","Chia1","Gkn1","Gkn2","Tff1","Tff2",
               "Muc5ac","Muc6","Mucl3","Psca","Lgals2","Cblif","Pdia2","Lypd8")
fib.genes <- c("Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a2","Col12a1",
               "Dcn","Bgn","Sparc","Lox","Grem1","Tnc","Hspg2","Pdgfra","Pdgfrb")
t.genes   <- c("Trac","Trbc1","Trbc2","Trdc","Tcrg-V6",
               "Cd3d","Cd3e","Cd3g","Cd247","Lef1","Zap70","Itk","Skap1","Rorc")
b.genes   <- c("Ighm","Ighd","Igkc","Iglc3","Cd79a","Cd79b","Ms4a1","Mzb1","Spib","Fcrla")

blacklist <- unique(c(epi.genes, fib.genes, t.genes, b.genes))

all_genes   <- rownames(seu)
m_gene_pool <- setdiff(all_genes, blacklist)

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  features = m_gene_pool
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_wo_blacklist.csv"))

# remove contamination, then re-cluster
seu <- subset(seu, subset = seurat_clusters %in% c(15, 20), invert = TRUE) # contamination
seu <- recluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))

# find markers, again, excluding obvious non-mye genes 
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  features = m_gene_pool
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_wo_blacklist_2.csv"))

# Check markers by feature plots
features <- readLines("aux_data/gene_set/additional/33_myeloid_markers.txt")
sapply(features, save_fp, seu, fp_path)

add_feat <- "Il1b"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_unfil.png", path = plot_path, width = 15, height = 3, units = "in", dpi = 150)

# adjust resolution if needed
# seu <- FindClusters(seu, resolution = 2, verbose = FALSE)

# save dim plots
DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters")
ggsave("cluster.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 45, height = 45, units = "mm")

# Add celltype annotation and save the Seurat object
seu$celltype2 <- ""
seu$celltype2[seu$seurat_clusters %in% c(10)] <- "Cyc.-Mye."
seu$celltype2[seu$seurat_clusters %in% c(4,7)] <- "Mac."
seu$celltype2[seu$seurat_clusters %in% c(1,3)] <- "Mono"
seu$celltype2[seu$seurat_clusters %in% c(5)] <- "TAM1"
seu$celltype2[seu$seurat_clusters %in% c(2,8)] <- "TAM2"
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "cDC2"
seu$celltype2[seu$seurat_clusters %in% c(12)] <- "Mig.-DC"
seu$celltype2[seu$seurat_clusters %in% c(14)] <- "moDC"
seu$celltype2[seu$seurat_clusters %in% c(15)] <- "pDC"
seu$celltype2[seu$seurat_clusters %in% c(0,9,11)] <- "TAN"
seu$celltype2[seu$seurat_clusters %in% c(13)] <- "Mast"
seu$celltype2 <- factor(seu$celltype2, levels = c("Cyc.-Mye.", "Mac.", "Mono", "TAM1", "TAM2",
                                                  "cDC2", "Mig.-DC", "moDC", "pDC", "TAN", "Mast"))
DimPlot(seu, group.by = "celltype2", cols = pal$myeloid, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, width = 35, height = 30, units = "mm")

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  features = m_gene_pool
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

features <- c(
  "Top2a", "Ankle1",
  "Pla2g1b", "Stab1",
  "Ccr2", "Vcan",
  "Il1a", "Cxcl2",
  "Apoe", "Mrc1",
  "H2-Ab1", "Cd74",
  "Ccr7", "Fscn1",
  "Ms4a4c", "H2-DMa",
  "Siglech", "Cacna1e",
  "S100a9", "Retnlg",
  "Tpsab1", "Cpa3"
)
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 85, height = 45, units = "mm")

saveRDS(seu, file = file.path("RDSfiles", "seu_033_myeloid.RDS"))

# Cell numbers and proportion plot ----
seu <- readRDS("RDSfiles/seu_033_myeloid.RDS")
fraction <- "myeloid"

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
  scale_fill_manual(values = pal$myeloid) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 40, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$myeloid) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 35, units = "mm")

# EDA ----
add_feat <- "Tnfrsf23"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend()
ggsave(paste0(add_feat, ".pdf"), path = fp_path, width = 25, height = 30, units = "mm")

features <- c("Adgre1", "Itgam", "S100a9", "Ly6g")
sapply(features, save_fp, seu, fp_path)
