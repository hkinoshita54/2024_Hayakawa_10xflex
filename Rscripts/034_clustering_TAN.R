# Continued from 033_clustering_mye.R
# Without integration
analysis_step <- "034_clustering_TAN"

# load packages ----
library(tidyverse)
library(Seurat)
library(msigdbr)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(scico)

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

# Load files ----
pal <- readRDS("RDSfiles/color_palette.RDS")
seu_mye <- readRDS("RDSfiles/seu_033_myeloid.RDS")
seu <- subset(seu_mye, subset = celltype2 == "TAN")

# Clustering ----
## regress out chrY related genes
y_genes <- c(
  "Ddx3y","Eif2s3y","Kdm5d","Uty","Zfy1","Zfy2","Sry",
  "Rps4y2","Usp9y","Uba1y","Jarid1d","Smcy","Tspy"
)
y_genes <- intersect(y_genes, rownames(seu))

seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL

seu <- AddModuleScore_UCell(seu, features = list(chrY = y_genes))

seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)

hvg <- VariableFeatures(seu)
VariableFeatures(seu) <- setdiff(hvg, y_genes)

seu <- ScaleData(seu, features = VariableFeatures(seu), vars.to.regress = c("chrY_UCell"))
seu <- RunPCA(seu, npcs = 30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(seu, dims = 1:30)

DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster.png", path = plot_path, width = 3, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) & NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

## adjust resolution if needed
# seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

# find markers with all genes ----
## remove b cell related genes (and chrY genes) from DE analysis
all_genes   <- rownames(seu)
b.genes   <- c("Ighm","Ighd","Igkc","Iglc3","Cd79a","Cd79b","Ms4a1","Mzb1","Spib","Fcrla")
blacklist <- c(y_genes, b.genes)
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
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

# MDSC signature ----
mdsc_sig <- openxlsx2::read_xlsx("aux_data/aay6017_table_s5.xlsx") %>% 
  pull(x)
pmn <- c("Ngp","Ltf","Camp","Csf3r","Cxcr2","Retnlg","S100a8","S100a9")
ifn <- c("Rtp4","Rsad2","Oas3","Isg15","Ifitm2")
supp <- c("Arg2","Il4ra","Cd84","Fgl2","Entpd1","Socs3","Osm")
mdsc_sig <- list(
  MDSC = mdsc_sig,
  PMN = pmn,
  IFN = ifn,
  SUPP = supp
)
seu <- AddModuleScore_UCell(
  seu,
  features = mdsc_sig,
  chunk.size = 1000,
  ncores = 8,
)

# Feature plot for MDSC
features <- grep("_UCell$", colnames(seu[[]]), value = TRUE)
sapply(features, save_fp, seu, fp_path)

add_feat <- "MDSC_UCell"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

#Vln
VlnPlot(seu, "MDSC_UCell", group.by = "seurat_clusters") & NoLegend()
ggsave("mdsc_vln_cluster.png", path = plot_path, width = 3, height = 3, units = "in", dpi = 150)
VlnPlot(seu, "IFN_UCell", group.by = "seurat_clusters") & NoLegend()
ggsave("ifn_vln_cluster.png", path = plot_path, width = 3, height = 3, units = "in", dpi = 150)

# HALLMARK
H <- msigdbr(species    = "Mus musculus", category   = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)
seu$MDSC_UCell <- NULL
seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

## Heatmap
meta <- seu[[]]
score_cols <- grep("_UCell$", colnames(meta), value = TRUE)
avg_scores <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
mat <- t(scale(t(mat)))
colnames(mat) <- avg_scores$seurat_clusters

rng <- range(mat, na.rm = TRUE)
seq_cols <- scico(9, palette  = "vik")   
col_fun <- colorRamp2(
  breaks = seq(rng[1], rng[2], length.out = length(seq_cols)),
  colors = seq_cols
)

ht <- Heatmap(
  mat,
  name = "UCell_z_score",
  row_labels = gsub("_UCell", "", score_cols),
  cluster_columns = FALSE,
  col = col_fun
)

pdf(file.path(plot_path,"UCell_Hallmark_heatmap.pdf"), width = 5, height = 10.5)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

add_feat <- "Selplg"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

