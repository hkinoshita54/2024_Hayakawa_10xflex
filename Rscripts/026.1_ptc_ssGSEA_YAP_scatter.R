# Continued from 026_ptc_ssGSEA.R
# 2026-02-23
# correlation of YAP and other scores
analysis_step <- "026_ptc_ssGSEA"

# load packages ----
library(tidyverse)
library(Seurat)
library(msigdbr)
library(UCell)
library(ggpubr)
source("Rscripts/theme_panel.R")

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

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

# FeatureScatter with YAP ----
## TME_score
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_TME_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## EMT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_EMT.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_MET_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## HYPOXIA
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_HYPOXIA.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_HYPOXIA_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

## INFL_RESP
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_INFL_RESP.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_INFL_RESP_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with NFkB ----
## TME_score
FeatureScatter(seu, "TNFA_SIGNALING_VIA_NFKB_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_NFkB_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "TNFA_SIGNALING_VIA_NFKB_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_NFkB_TME_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with G2M ----
## TME_score
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_G2M_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_G2M_TME_no_p.pdf", path = plot_path, width = 30, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with nFeature_RNA ----
## TME_score
FeatureScatter(seu, "nFeature_RNA", "TME_score", group.by = "celltype3.1", cols = ptc_cols, pt.size = 0.2) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_nFeature_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)
