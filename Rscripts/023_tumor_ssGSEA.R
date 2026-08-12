# Continued from 021_clustering_epi_PC&PTC.R
# Without integration
analysis_step <- "023_tumor_ssGSEA"

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
pal <- readRDS("RDSfiles/color_palette.RDS")
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
Idents(seu) <- "celltype3"
DimPlot(seu, group.by = "celltype3", cols = pal$tum, label = TRUE, repel = TRUE) & NoAxes()

# UCell scoring ----
H <- msigdbr(species    = "Mus musculus", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

# C6 <- msigdbr(species    = "Mus musculus", collection = "C6")
# C6 <- split(C6$gene_symbol, C6$gs_name)

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
  group_by(celltype3) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype3
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

## plot only selected gene sets
rownames(mat)
gene_sets_sel <- c(
  "E2F_TARGETS_UCell",
  "G2M_CHECKPOINT_UCell",
  "ANGIOGENESIS_UCell",
  "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell",
  "HYPOXIA_UCell",
  "INFLAMMATORY_RESPONSE_UCell",
  "TNFA_SIGNALING_VIA_NFKB_UCell",
  "GLYCOLYSIS_UCell",
  "OXIDATIVE_PHOSPHORYLATION",
  "FATTY_ACID_METABOLISM"
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
  name = "UCell_z_score",
  row_labels = gsub("_UCell", "", rownames(mat_sel)),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun
)

pdf(file.path(plot_path,"UCell_Hallmark_heatmap_sel.pdf"), width = 5.5, height = 4)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()


# UCell scoring without Trans. and Ca. 4 ----
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
Idents(seu) <- "celltype3"
seu <- subset(seu, subset = celltype3 %in% c("Trans.", "Ca.4"), invert = TRUE)
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
  group_by(celltype3) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))
mat <- t(as.matrix(avg_scores[, score_cols]))
colnames(mat) <- avg_scores$celltype3
openxlsx2::write_xlsx(mat, row_names = TRUE, file = file.path(res_path, "H_UCell_2.xlsx"))
mat <- t(scale(t(mat)))

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

pdf(file.path(plot_path,"UCell_Hallmark_heatmap_2.pdf"), width = 6, height = 10.5)  # adjust size
draw(
  ht,
  heatmap_legend_side = "top",
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()


# output for scenic ----
out_dir <- file.path("out_data", analysis_step)
fs::dir_create(out_dir)

## load data
# seu <- readRDS("RDSfiles/seu_021_tumor.RDS")

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
