# Continued from 020_clustering_epi_2.R
# Without integration
analysis_step <- "022_epi_ssGSEA"

# load packages ----
library(tidyverse)
library(Seurat)
library(msigdbr)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(scico)
library(grid)
source("Rscripts/theme_panel.R")

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# load file and check dimplot ----
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))
Idents(seu) <- "celltype2"
DimPlot(seu, group.by = "celltype2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
seu_cp <- seu

H <- msigdbr(species    = "Mus musculus", category   = "H")
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
  group_by(celltype2) %>%
  dplyr::summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype2
mat <- t(scale(t(mat)))
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "H_UCell.xlsx"))

# rng <- range(mat, na.rm = TRUE)
seq_cols <- scico(9, palette  = "vik")  
col_fun <- colorRamp2(
  breaks = seq(-2, 2, length.out = length(seq_cols)),
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
# lgd <- Legend(
#   title = "Z-score",
#   col_fun = col_fun,
#   title_gp = gpar(fontsize = 7),
#   labels_gp = gpar(fontsize = 6)
# )
pdf(file.path(plot_path,"UCell_Hallmark_heatmap.pdf"), width = 4.5, height = 5.5)  # adjust size
draw(
  ht,
  # heatmap_legend_list = list(lgd),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()

# UCell scoring yap signatures ----
## yap signature
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
  group_by(celltype2) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype2
mat <- t(scale(t(mat)))
# openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "YAP_UCell.xlsx"))

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

pdf(file.path(plot_path,"UCell_yap_signatures_heatmap.pdf"), width = 5, height = 2.6)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

# UCell scoring with Stem, Prog, Epi-PC, Epi-PTC ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")
Idents(seu) <- "celltype2"
seu <- subset(seu, subset = celltype2 %in% c("Stem", "Prog.", "Epi-PC", "Epi-PTC"))
seu_cp <- seu

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
  group_by(celltype2) %>%
  dplyr::summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype2
mat <- t(scale(t(mat)))
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "H_UCell_2.xlsx"))

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
  cluster_rows = TRUE,
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
pdf(file.path(plot_path,"UCell_Hallmark_heatmap_4cols.pdf"), width = 3.6, height = 5.5)  # adjust size
draw(
  ht,
  heatmap_legend_list = list(lgd),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()

## plot only selected gene sets
rownames(mat)
gene_sets_sel <- c(
  # "ANGIOGENESIS_UCell",
  "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell",
  "HYPOXIA_UCell",
  "INFLAMMATORY_RESPONSE_UCell",
  "TNFA_SIGNALING_VIA_NFKB_UCell",
  # "GLYCOLYSIS_UCell",
  "OXIDATIVE_PHOSPHORYLATION",
  "FATTY_ACID_METABOLISM",
  "E2F_TARGETS_UCell",
  "G2M_CHECKPOINT_UCell"
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
pdf(file.path(plot_path,"UCell_Hallmark_heatmap_sel_2.pdf"), width = 3.1, height = 1.5)  # adjust size
draw(
  ht,
  heatmap_legend_list = list(lgd),
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
  group_by(celltype2) %>%
  dplyr::summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
mat <- t(scale(t(mat)))
colnames(mat) <- avg_scores$celltype2
# openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "YAP_UCell_2.xlsx"))

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
pdf(file.path(plot_path,"UCell_yap_signatures_4clos.pdf"), width = 1.5, height = 0.85)  # adjust size
draw(
  ht,
  # heatmap_legend_list = list(lgd),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()
