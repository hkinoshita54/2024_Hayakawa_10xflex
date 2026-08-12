# Continued from 026_ptc_ssGSEA.R
# 2026-06-06
# compute scores of metaprograms from human GC scRNA-seq

# Make directories ----
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"
analysis_step <- "026.2_ptc_metaprogram"
plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

# load packages ----
library(tidyverse)
library(Seurat)
library(msigdbr)
library(UCell)
library(ggpubr)
source("Rscripts/theme_panel.R")
source(file.path(ws, "common_Rscripts/helpers.R"))

# load file and check dimplot ----
ptc_cols <- c(PTC1 = "#CC79A7", PTC2 = "#7F7F7F", PTC3 = "#E69F00", PTC4 = "#009E73")
seu <- readRDS("RDSfiles/seu_025_ptc.RDS")
Idents(seu) <- "celltype3.1"
DimPlot(seu, group.by = "celltype3.1", cols = ptc_cols) & NoAxes()

## load gene sets
H <- msigdbr(species    = "Mus musculus", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list.RDS")
yap_list <- yap_list[4:5]

MP_list <- readRDS("RDSfiles/MP_list_mm.RDS")

# UCell scoring ----
seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = MP_list,
  chunk.size = 1000,
  ncores = 8,
)

## make ucell scores into z-scores

### change names for brevity
colnames(seu@meta.data)[colnames(seu@meta.data) == "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"] <- "EMT_UCell"
colnames(seu@meta.data)[colnames(seu@meta.data) == "INFLAMMATORY_RESPONSE_UCell"] <- "INFL_RESP_UCell"

ucell_cols <- grep("_UCell$", colnames(seu[[]]), value = TRUE)

for (cc in ucell_cols) {
  seu[[cc]] <- as.numeric(scale(seu@meta.data[[cc]]))
}

## create a composite score
seu[["TME_score"]] <- rowMeans(
  seu@meta.data[, c(
    "EMT_UCell",
    "HYPOXIA_UCell",
    "INFL_RESP_UCell"
  )],
  na.rm = TRUE
)

sapply(ucell_cols, save_fp, seu, fp_path, pt.size = 6)
sapply(MP_list$core9, save_fp, seu, fp_path, pt.size = 6)

features <- c("Lamc2", "Lamc3", "Rhoc", "Tnfrsf12a", "Krt7", "S100a14", "Sfn")
sapply(features, save_fp, seu, fp_path, pt.size = 6)

# FeatureScatter with YAP ----

## MP4 and TME_score
FeatureScatter(seu, "MP_4_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_TME_score.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_4_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_TME_score_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP4 and G2M
FeatureScatter(seu, "MP_4_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_G2M.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_4_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_G2M_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP4 and YAP
FeatureScatter(seu, "MP_4_UCell", "YAP_TARGETS_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_YAP.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_4_UCell", "YAP_TARGETS_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP4_YAP_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP_5 and TME_score
FeatureScatter(seu, "MP_5_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP5_TME_score.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_5_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP5_TME_score_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP_5 and G2M
FeatureScatter(seu, "MP_5_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP5_G2M.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_5_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP5_G2M_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP14 and TME_score
FeatureScatter(seu, "MP_14_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP14_TME_score.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_14_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP14_TME_score_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## MP14 and G2M
FeatureScatter(seu, "MP_14_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP14_G2M.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "MP_14_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_MP14_G2M_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## core9 and TME_score
FeatureScatter(seu, "core9_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_core9_TME_score.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "core9_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_core9_TME_score_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## core9 and G2M
FeatureScatter(seu, "core9_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_core9_G2M.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "core9_UCell", "G2M_CHECKPOINT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_core9_G2M_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)


# dotplot for grant ----
## 2026-07-08
seu$invasion <- case_when(
  seu$celltype3.1 == "PTC3" ~ "Invasive",
  seu$celltype3.1 %in% c("PTC1", "PTC2") ~ "Non_inv",
  TRUE ~ NA_character_,
)
seu <- subset(seu, subset = celltype3.1 == "PTC4", invert = TRUE)
seu$invasion <- factor(seu$invasion, levels = c("Non_inv", "Invasive"))
features <- c(
  "Areg",
  "Ets2",
  "Lamc2",
  "Plaur",
  "Anxa1",
  "Krt7"
)
DotPlot(seu, group.by = "invasion", features = features, dot.scale = 2.3) + 
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
  )+
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) +
  coord_flip()
ggsave("dotplot_inv_mp.pdf", path = plot_path, width = 35, height = 32, units = "mm")
