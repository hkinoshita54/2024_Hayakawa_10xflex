# doing this AFTER 052
analysis_step <- "064_reference_ingest_PTC"

# load packages ----
library(tidyverse)
library(Seurat)
source("Rscripts/theme_panel.R")

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
tum <- readRDS("RDSfiles/seu_025_ptc.RDS")
table(seu_all$orig.ident, seu_all$celltype1)
levels(seu_all$celltype2)
levels(tum$celltype3.1)

# subset to PTC, with the exception for Squam., Bcell, Glial and Neural ----
## drop Stem, Pit, EEC, Tuft, Adipo., Meso. which we don't analyze later
seu <- subset(seu_all, 
              subset = orig.ident %in% c("PTC") |
                celltype1 %in% c("Squam.", "Bcell", "Adipo.", "Glial", "Neural", "Meso.")
              )
seu <- subset(seu, subset = celltype1 %in% c("Stem", "Pit", "EEC", "Tuft"), invert = TRUE)
seu$celltype1 <- droplevels(seu$celltype1)
seu$celltype2 <- droplevels(seu$celltype2)
table(seu$cellgroup)
table(seu$celltype1)
table(seu$celltype2)

# Group annotations ----
## rename or group
seu$level1 <- case_when(
  seu$celltype1 %in% c("Epi-PC&PTC") ~ "Tumor",
  TRUE ~ as.character(seu$celltype1)
)
seu$level1 <- factor(seu$level1, levels = c("Squam.","Tumor","Tcell","Bcell","Mye.","Fibro.","EC","Myo.","Peri.","Adipo.","Glial","Neural","Meso."))
table(seu$level1)

seu$level2 <- as.character(seu$celltype2)

all(names(seu$level2[seu$level2=="Epi-PTC"]) == names(tum$celltype3.1))
seu$level2[seu$level2=="Epi-PTC"] <- as.character(tum$celltype3.1)
levels <- levels(seu$celltype2)
levels <- c("Squam.", levels(tum$celltype3.1), levels[3:42])
seu$level2 <- factor(seu$level2, levels = levels)
table(seu$level2)

# downsample b cells from PC so PC B cells not dominate the reference ----
set.seed(2026)

## get cell names
b_cells <- WhichCells(seu, expression = celltype1 == "Bcell")
b_pc   <- WhichCells(seu, expression = celltype1 == "Bcell" & orig.ident == "PC")
b_nonpc <- WhichCells(seu, expression = celltype1 == "Bcell" & orig.ident != "PC")

## downsample to 600
b_pc_keep <- sample(b_pc, 600) 
b_keep <- c(b_pc_keep, b_nonpc)
cells_keep <- c(setdiff(Cells(seu), b_cells), b_keep)

## subset Seurat object
seu <- subset(seu, cells = cells_keep)
table(seu$orig.ident, seu$celltype1)

# make the object clean and save ----
table(seu$cellgroup)
table(seu$level1)
table(seu$level2)

saveRDS(seu, "RDSfiles/seu_064_ref_ptc.RDS")

# plots for record ----
seu <- cluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 45, height = 45, units = "mm")

DimPlot(seu, group.by = "level1", cols = "glasbey", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes()
ggsave("level1.pdf", path = plot_path, width = 60, height = 45, units = "mm")

DimPlot(seu, group.by = "level2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes()
ggsave("level2.pdf", path = plot_path, width = 75, height = 45, units = "mm")



# export data for python ----
out_dir <- file.path("out_data", analysis_step)
fs::dir_create(out_dir)

# save metadata table
seu$barcode <- colnames(seu)
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
meta <- seu@meta.data %>% select(barcode, UMAP_1, UMAP_2, orig.ident, cellgroup, level1, level2)
write.csv(meta, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

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

