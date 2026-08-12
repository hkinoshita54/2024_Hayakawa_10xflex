# Continued from 050_imm_combined.R
analysis_step <- "050.3_imm_ssGSEA_2"

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

# Helper functions ----
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

pal <- readRDS("RDSfiles/color_palette.RDS")

# load file and check dimplot ----
seu <- readRDS("RDSfiles/seu_050_imm_combined.RDS")
Idents(seu) <- "celltype2"
DimPlot(seu, group.by = "celltype2", label = TRUE, repel = TRUE,
        cols = c(pal$T_cells, pal$B_cells, pal$myeloid)) & NoAxes()
seu_cp <- seu

# UCell scoring of MDSC signatures ----
# MDSC ----
mdsc_sig <- openxlsx2::read_xlsx("aux_data/aay6017_table_s5.xlsx") %>% 
  pull(x)
mdsc_sig <- list(
  MDSC = mdsc_sig
)
seu <- AddModuleScore_UCell(
  seu,
  features = mdsc_sig,
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
mat <- t(scale(t(mat)))
colnames(mat) <- avg_scores$celltype2

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

pdf(file.path(plot_path,"UCell_MDSC_sig_heatmap.pdf"), width = 8, height = 2.8)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

# feature plot ----
sapply(score_cols, save_fp, seu, fp_path)

add_feat <- "Msr1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# dotplot ----
DotPlot(seu, group.by = "celltype2", features = "MDSC_UCell") & RotatedAxis()
ggsave("mdsc_dotplot.png", path = plot_path, width = 4, height = 8, units = "in", dpi = 150)

## export mdsc_sig to a table
max_len <- max(lengths(mdsc_sig))
df <- as.data.frame(lapply(mdsc_sig, function(x) {
  c(x, rep(NA, max_len - length(x)))
}))
df <- t(df)
openxlsx2::write_xlsx(
  df, 
  col_names = FALSE,
  row_names = TRUE,
  na.strings = "",
  file = file.path(res_path, "mdsc_sig.xlsx"),
  
)
