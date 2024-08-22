# Continued from 020_clustering_epi.R and pySCENIC output (nig-sc, 5th iteration)
analysis_step <- "200_epi_scenic"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# load data ----
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))

# auc_mtx is exported from pySCENIC (5th iteration)
auc_mtx <- read_tsv(file = file.path("int_data", "auc_mtx.txt")) %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
seu[["scenicAUC"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "scenicAUC"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)

DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("regulon_celltype.xlsx")))

Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
markers <- FindAllMarkers(seu, only.pos = TRUE)

add_feat <- "Mki67"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = TRUE, repel = TRUE
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)
rm(add_feat)

Idents(seu) <- "celltype"
saveRDS(seu, file = file.path("RDSfiles", "seu_200_epi_scenic.RDS"))
