# Continued from 050_imm_combined.R
analysis_step <- "050.2_imm_ssGSEA"

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
load("RDSfiles/immune_signatures.RData")
seu <- readRDS("RDSfiles/seu_050_imm_combined.RDS")
seu <- cluster(seu, npcs = 50, res = 1)
Idents(seu) <- "celltype2"
DimPlot(seu, group.by = "celltype2", label = TRUE, repel = TRUE,
        cols = c(pal$T_cells, pal$B_cells, pal$myeloid)) & NoAxes()
seu_cp <- seu

#####
# this is test of scGSVA with small sample

# UCell scoring of HALLMARKS ----
H <- msigdbr(species    = "Mus musculus", category   = "H")
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
  col = col_fun
)

pdf(file.path(plot_path,"UCell_Hallmark_heatmap.pdf"), width = 8, height = 10.5)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()


# UCell scoring of immune_signatures ----
seu <- seu_cp
seu <- AddModuleScore_UCell(
  seu,
  features = immune_signatures,
  chunk.size = 1000,
  ncores = 8,
)

# Heatmap, zscores by row
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

pdf(file.path(plot_path,"UCell_imm_sig_heatmap.pdf"), width = 9, height = 7)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

# Heatmap, zscores by column
meta <- seu[[]]
score_cols <- grep("_UCell$", colnames(meta), value = TRUE)
avg_scores <- meta %>%
  group_by(celltype2) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
mat <- scale(mat)
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

pdf(file.path(plot_path,"UCell_imm_sig_heatmap_z_by_col.pdf"), width = 9, height = 7)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()


# Cell numbers and proportion plot ----
fraction <- "imm_combined"

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
  scale_fill_manual(values = c(pal$T_cells, pal$B_cells, pal$myeloid)) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 3.6, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c(pal$T_cells, pal$B_cells, pal$myeloid)) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 3.6, height = 3.2, units = "in", dpi = 150)


## export immune_signatures to a table
max_len <- max(lengths(immune_signatures))
df <- as.data.frame(lapply(immune_signatures, function(x) {
  c(x, rep(NA, max_len - length(x)))
}))
df <- df[,20:1]
df <- t(df)
openxlsx2::write_xlsx(
  df, 
  col_names = FALSE,
  row_names = TRUE,
  na.strings = "",
  file = file.path(res_path, "immune_signatures.xlsx"),
  
)
