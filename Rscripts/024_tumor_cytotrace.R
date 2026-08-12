# Continued from 021_clustering_epi_PC_PTC.R
# Without integration
analysis_step <- "024_tumor_cytotrace"

# load packages ----
library(tidyverse)
library(Seurat)
library(CytoTRACE)
library(CytoTRACE2)
# library(monocle3)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# load data
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

# CytoTRACE ----
Idents(seu) <- "celltype3"
counts_matrix <- LayerData(seu, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu@meta.data$celltype3)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu@meta.data$celltype3)
names(pheno) <- colnames(seu)
umap_df <- seu@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = pheno, emb = umap_df, outputDir = res_path)

# CytoTRACE2 ----
seu <- cytotrace2(seu, is_seurat = TRUE, slot_type = "counts", species = 'mouse')
add_feat <- "CytoTRACE2_Relative"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = FALSE, repel = TRUE
) + scale_color_viridis_c(option = "B") + NoAxes()
ggsave(paste0(add_feat, ".png"), path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

VlnPlot(seu, features = add_feat, cols = pal$tum, pt.size = 0)
ggsave("vln_cytotrace2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
