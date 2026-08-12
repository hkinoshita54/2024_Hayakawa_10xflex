# Continued from 025_clustering_epi-PTC.R
# 2026-01-10
analysis_step <- "026_ptc_ssGSEA"

# load packages ----
library(tidyverse)
library(Seurat)
library(msigdbr)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(scico)
# library(CytoTRACE)
library(CytoTRACE2)
library(grid)
source("Rscripts/theme_panel.R")

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# load file and check dimplot ----
ptc_cols <- c(PTC1 = "#CC79A7", PTC2 = "#7F7F7F", PTC3 = "#E69F00", PTC4 = "#009E73")
seu <- readRDS("RDSfiles/seu_025_ptc.RDS")
seu_cp <- seu
Idents(seu) <- "celltype3.1"
DimPlot(seu, group.by = "celltype3.1", cols = ptc_cols) & NoAxes()

H <- msigdbr(species    = "Mus musculus", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list.RDS")
yap_list <- yap_list[4:5]

# UCell scoring ----
seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

# Heatmap
meta <- seu[[]]
score_cols <- grep("_UCell$", colnames(meta), value = TRUE)
avg_scores <- meta %>%
  group_by(celltype3.1) %>%
  dplyr::summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype3.1
mat <- t(scale(t(mat)))
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "H_UCell.xlsx"))


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
  cluster_rows = FALSE,
  col = col_fun
)

pdf(file.path(plot_path,"UCell_Hallmark_heatmap.pdf"), width = 5, height = 10.5)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

## plot only selected gene sets
rownames(mat)
gene_sets_sel <- c(
  "E2F_TARGETS_UCell",
  "G2M_CHECKPOINT_UCell",
  "FATTY_ACID_METABOLISM",
  "OXIDATIVE_PHOSPHORYLATION",
  # "PROTEIN_SECRETION",
  "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell",
  # "TGF_BETA_SIGNALING",
  "APICAL_JUNCTION",
  "HYPOXIA_UCell",
  "INFLAMMATORY_RESPONSE"
)
mat_sel <- as.data.frame(mat)[gene_sets_sel,] %>% as.matrix

rng <- range(mat_sel, na.rm = TRUE)
seq_cols <- scico(9, palette  = "vik") 
col_fun <- colorRamp2(
  breaks = seq(rng[1], rng[2], length.out = length(seq_cols)),
  colors = seq_cols
)

ht <- Heatmap(
  mat_sel,
  row_labels = gsub("_UCell", "", rownames(mat_sel)),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  show_heatmap_legend = FALSE
)
lgd <- Legend(
  title = "Z-score",
  col_fun = col_fun,
  title_gp = gpar(fontsize = 7),
  labels_gp = gpar(fontsize = 6)
)
pdf(file.path(plot_path,"UCell_Hallmark_heatmap_sel_2.pdf"), width = 2.8, height = 1.3)  # adjust size
draw(
  ht,
  # heatmap_legend_list = list(lgd),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()

# UCell scoring yap signatures ----
seu <- seu_cp
seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

# Heatmap
meta <- seu[[]]
score_cols <- grep("_UCell$", colnames(meta), value = TRUE)
avg_scores <- meta %>%
  group_by(celltype3.1) %>%
  dplyr::summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
mat <- t(scale(t(mat)))
colnames(mat) <- avg_scores$celltype3.1
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "YAP_UCell_2.xlsx"))

rng <- range(mat, na.rm = TRUE)
seq_cols <- scico(9, palette  = "vik")  
col_fun <- colorRamp2(
  breaks = seq(rng[1], rng[2], length.out = length(seq_cols)),
  colors = seq_cols
)

ht <- Heatmap(
  mat,
  row_labels = gsub("_UCell", "", rownames(mat)),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  show_heatmap_legend = FALSE
)
lgd <- Legend(
  title = "Z-score",
  col_fun = col_fun,
  title_gp = gpar(fontsize = 7),
  labels_gp = gpar(fontsize = 6)
)
pdf(file.path(plot_path,"UCell_yap_signatures.pdf"), width = 1.5, height = 0.7)  # adjust size
draw(
  ht,
  # heatmap_legend_list = list(lgd),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()

# CytoTRACE ----
# counts_matrix <- LayerData(seu, assay='RNA', layer='counts') %>% as.data.frame()
# obj_cell_type_anno <- as.data.frame(seu@meta.data$celltype3)
# results <- CytoTRACE(counts_matrix, ncores = 4)
# pheno <- as.character(seu@meta.data$celltype3.1)
# names(pheno) <- colnames(seu)
# umap_df <- seu@reductions$umap@cell.embeddings
# plotCytoTRACE(results, phenotype = pheno, emb = umap_df, outputDir = res_path)

# CytoTRACE2 ----
counts <- LayerData(seu, assay='RNA', layer='counts')
ct2 <- cytotrace2(counts, species = "mouse")
seu$CytoTRACE2_Relative <- ct2$CytoTRACE2_Relative
add_feat <- "CytoTRACE2_Relative"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = FALSE, repel = TRUE
) + scale_color_viridis_c(option = "B") + NoAxes()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 3, units = "in", dpi = 150)

VlnPlot(seu, features = add_feat, cols = ptc_cols, pt.size = 0)  &
  theme_panel() & NoLegend() &
  labs(title = "CytoTRACE2", x = NULL, y = "relative_score") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 30, height = 40, units = "mm")


# output for scenic ----
out_dir <- file.path("out_data", analysis_step)
fs::dir_create(out_dir)

## load data
# seu <- readRDS("RDSfiles/seu_025_ptc.RDS")

## convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

## save metadata table
seu$barcode <- colnames(seu)
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
seu <- FindVariableFeatures(seu, nfeatures = 8000)
hvgs <- VariableFeatures(seu)
tfs <- readLines("pyscenic_docker/aux_data/allTFs_mm.txt")
tfs <- tfs[tfs %in% rownames(seu)]
genes_for_scenic <- union(hvgs, tfs)
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
counts_matrix <- counts_matrix[genes_for_scenic,]
writeMM(counts_matrix, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)
