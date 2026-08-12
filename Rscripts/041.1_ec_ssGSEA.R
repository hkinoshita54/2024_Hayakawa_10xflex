# Continued from 041_clustering_ec.R
# Without integration
analysis_step <- "041.1_ec_ssGSEA"

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

# load file and check dimplot ----
source("Rscripts/theme_panel.R")
pal <- readRDS("RDSfiles/color_palette.RDS")
seu <- readRDS("RDSfiles/seu_041_ec.RDS")
Idents(seu) <- "celltype2"
DimPlot(seu, group.by = "celltype2", cols = pal$endothelial, label = TRUE, repel = TRUE) & NoAxes()
seu_cp <- seu

# UCell scoring ----
H <- msigdbr(species    = "Mus musculus", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

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
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype2
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
  col = col_fun
)

pdf(file.path(plot_path,"UCell_Hallmark_heatmap.pdf"), width = 7, height = 10.5)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

# UCell scoring adhesion and chemokines ----
adhesion <- c("Icam1", "Icam2", "F11r", "Selp", "Sele")
chemokines <- c("Csf3", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5")
secreted <- c("Lrg1",  "Csf3",  "Pdgfa",  "Pdgfb", "Sema3f","Grn", "Tgfb1")
ec_scores <- list(
  adhesion = adhesion,
  chemokines = chemokines,
  secreted = secreted
)

seu <- seu_cp
Idents(seu) <- "celltype2"
seu <- AddModuleScore_UCell(
  seu,
  features = ec_scores,
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
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "H_UCell_2.xlsx"))
mat <- t(scale(t(mat)))

rng <- range(mat, na.rm = TRUE)
seq_cols <- viridis(9, option = "D")  
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

pdf(file.path(plot_path,"UCell_ec_scores_heatmap.pdf"), width = 6, height = 3)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

VlnPlot(seu, features = grep("_UCell$", colnames(seu[[]]), value = TRUE), 
        cols = pal$endothelial,
        ncol = 3, pt.size = 0)
ggsave("ec_signatures.png", path = plot_path, width = 8, height = 3, units = "in", dpi = 150)

signature = "secreted"
VlnPlot(seu, features = paste0(signature, "_UCell"), cols = pal$endothelial, pt.size = 0)  &
  theme_panel() & NoLegend() &
  labs(title = signature, x = NULL, y = "UCell_score") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", signature, ".pdf"), path = plot_path, width = 45, height = 45, units = "mm")
