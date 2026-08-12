# doing this AFTER 052
analysis_step <- "061_reference_spatial_WT"

# load packages ----
library(tidyverse)
library(Seurat)

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
seu_all <- readRDS("RDSfiles/seu_052_all_combined.RDS")
# tum <- readRDS("RDSfiles/seu_021_tumor.RDS")
table(seu_all$orig.ident, seu_all$celltype1)
levels(seu_all$celltype2)
# levels(tum$celltype3)

# subset to WT and PT, with the exception for Squam., Bcell, Glial and Neural ----
## drop Tuft, Pit-PT, Adipo., Meso. which we don't analyze later
seu <- subset(seu_all, subset = orig.ident %in% c("WT", "PT") | celltype1 %in% c("Squam.", "Bcell", "Glial", "Neural"))
seu <- subset(seu, subset = celltype1 %in% c("Tuft", "Pit-PT", "Adipo.", "Meso."), invert = TRUE)
seu$celltype1 <- droplevels(seu$celltype1)
seu$celltype2 <- droplevels(seu$celltype2)
table(seu$orig.ident, seu$celltype1)
table(seu$celltype2)

# Group annotations ----
## rename or group
seu$celltype4 <- case_when(
  seu$celltype1 %in% c("Myo.", "Peri.") ~ "Myo.",
  seu$celltype1 %in% c("Glial", "Neural") ~ "ENS",
  seu$celltype2 %in% c("TAN") ~ "TAN",
  seu$celltype2 %in% c("Cyc.-Mye.","Mac.","Mono","TAM1","TAM2","cDC2","Mig.-DC","moDC","pDC", "Mast") ~ "MNP",
  TRUE ~ as.character(seu$celltype1)
)
table(seu$celltype4)
seu <- subset(seu, subset = celltype4 == "TAN", invert = TRUE)

# downsample b cells from PC so PC B cells not dominate the reference ----
set.seed(2026)

## get cell names
b_cells <- WhichCells(seu, expression = celltype4 == "Bcell")
b_pc   <- WhichCells(seu, expression = celltype4 == "Bcell" & orig.ident == "PC")
b_nonpc <- WhichCells(seu, expression = celltype4 == "Bcell" & orig.ident != "PC")

## downsample to 600
b_pc_keep <- sample(b_pc, 600) 
b_keep <- c(b_pc_keep, b_nonpc)
cells_keep <- c(setdiff(Cells(seu), b_cells), b_keep)

## subset Seurat object
seu <- subset(seu, cells = cells_keep)
table(seu$orig.ident, seu$celltype4)

# make the object clean and save ----
seu$celltype4 <- factor(
  seu$celltype4, 
  levels = c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.","EEC", "Squam.", "Bcell", "Tcell", "MNP", "EC", "Fibro.", "Myo.", "ENS")
)

# plots for record ----
seu <- cluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) + NoAxes()
ggsave("sample.png", path = plot_path, width = 3.3, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype4", cols = "polychrome") & NoAxes()
ggsave("celltype4.png", path = plot_path, width = 3.4, height = 3, units = "in", dpi = 150)

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
cts <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(cts, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(cts)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)

# select genes for rctd ----
## 2000 hvg + top 50 marker genes for each cluster - mito, ribo, cell cycle
hvg <- VariableFeatures(seu)

Idents(seu) <- "celltype4"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers_use <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  pull(gene) %>% 
  unique()

genes_use <- unique(c(hvg, markers_use))
genes_use <- genes_use[!grepl("^mt-|^Rpl|^Rps", genes_use, ignore.case=FALSE)]
# load("RDSfiles/cellcycle_genes_mouse.Rdata")
# genes_use <- setdiff(genes_use, unique(c(s.genes, g2m.genes)))

## save for later use as a reference
metadata <- seu@meta.data %>% 
  select("orig.ident", "cellgroup", "celltype1", "celltype2", "celltype4")
cts <- LayerData(seu, assay = 'RNA', layer = 'counts')
seu_ref <- CreateSeuratObject(cts, meta.data = metadata)
saveRDS(seu_ref, "RDSfiles/seu_061_ref_WT.RDS")

seu_rctd <- seu_ref[genes_use,]
cts_rctd <- LayerData(seu_rctd, assay = 'RNA', layer = 'counts')
seu_rctd <- CreateSeuratObject(cts_rctd, meta.data = metadata)
saveRDS(seu_rctd, "RDSfiles/seu_061_ref_WT_rctd.RDS")
