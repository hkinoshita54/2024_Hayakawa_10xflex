# doing this AFTER 052
analysis_step <- "054_combine_for_deconv"

# load packages ----
library(tidyverse)
library(Seurat)
library(speckle)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

# load data ----
pal <- readRDS("RDSfiles/color_palette.RDS")
seu <- readRDS("RDSfiles/seu_052_all_combined.RDS")
levels(seu$celltype1)
levels(seu$celltype2)

# Group annotations ----
## group non-tumor epithelial cells, peri and myo, glial and neural
seu$celltype_deconv <- case_when(
  seu$celltype1 %in% c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft", "Pit-PT") ~ "NT-Epi",
  seu$celltype1 %in% c("Epi-PC&PTC") ~ "Tum-Epi",
  seu$celltype1 %in% c("Myo.", "Peri.") ~ "Myo.",
  seu$celltype1 %in% c("Glial", "Neural") ~ "ENS",
  seu$celltype2 %in% c("TAN") ~ "TAN",
  seu$celltype2 %in% c("Cyc.-Mye.","Mac.","Mono","TAM1","TAM2","cDC2","Mig.-DC","moDC","pDC", "Mast") ~ "MNP",
  TRUE ~ as.character(seu$celltype1)
)

# subset to WT and PTC, with the exception for Squam. and Bcell ----
## also drop adipo. and meso. which we don't analyze later
seu <- subset(seu, subset = orig.ident %in% c("WT", "PTC") | celltype_deconv %in% c("Squam.", "Bcell"))
seu <- subset(seu, subset = celltype_deconv %in% c("Adipo.", "Meso."), invert = TRUE)
table(seu$celltype_deconv)

# downsample b cells from PC so PC B cells not dominate the reference ----
set.seed(123)

## get cell names
b_cells <- WhichCells(seu, expression = celltype_deconv == "Bcell")
b_pc   <- WhichCells(seu, expression = celltype_deconv == "Bcell" & orig.ident == "PC")
b_nonpc <- WhichCells(seu, expression = celltype_deconv == "Bcell" & orig.ident != "PC")

## downsample to 600
b_pc_keep <- sample(b_pc, 600) 
b_keep <- c(b_pc_keep, b_nonpc)
cells_keep <- c(setdiff(Cells(seu), b_cells), b_keep)

## subset Seurat object
seu <- subset(seu, cells = cells_keep)
table(seu$orig.ident, seu$celltype_deconv)

# make the object clean and save ----
seu$celltype_deconv <- factor(
  seu$celltype_deconv, 
  levels = c("NT-Epi", "Tum-Epi", "Squam.", "Bcell", "Tcell", "MNP", "TAN", "EC", "Fibro.", "Myo.", "ENS")
)
seu
metadata <- seu@meta.data %>% 
  select("orig.ident", "cellgroup", "celltype1", "celltype2", "celltype_deconv")
cts <- LayerData(seu, assay = 'RNA', layer = 'counts')
seu_ref <- CreateSeuratObject(cts, meta.data = metadata)

table(seu_ref$celltype_deconv)
table(seu_ref$orig.ident, seu_ref$celltype_deconv)
saveRDS(seu, "RDSfiles/seu_054_for_deconv.RDS")

# plots for record ----
seu <- cluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) + NoAxes()
ggsave("sample.png", path = plot_path, width = 3.3, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype_deconv", cols = "polychrome") & NoAxes()
ggsave("celltype_deconv.png", path = plot_path, width = 3.4, height = 3, units = "in", dpi = 150)

# export data for python ----
out_dir <- file.path("out_data", analysis_step)
fs::dir_create(out_dir)

# save metadata table
seu$barcode <- colnames(seu)
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(counts_matrix, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)
