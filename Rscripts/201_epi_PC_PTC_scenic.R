# Continued from 020_clustering_epi.R and pySCENIC output (nig-sc, 5th iteration)
analysis_step <- "201_epi_PC_PTC_scenic"

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
seu <- readRDS(file.path("RDSfiles", "seu_021_epi_PC_PTC.RDS"))

# auc_mtx is exported from pySCENIC
# first check the auc_mtx from the analysis of whole epithelium
plot_path <- file.path("plot", analysis_step, "auc_mtx_from_Epithelial")
res_path <- file.path("result", analysis_step, "auc_mtx_from_Epithelial")
fs::dir_create(c(plot_path, res_path))

auc_mtx <- read_tsv(file = file.path("int_data", "pyscenic_nig-sc", "Epithelial", "auc_mtx.txt")) %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[, colnames(seu)]
seu[["scenicAUC"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "scenicAUC"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)

DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype_2", label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype_2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, "regulon_celltype.xlsx"))

Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
markers <- FindAllMarkers(seu, only.pos = TRUE)

DefaultAssay(seu) <- "scenicAUC"
add_feat <- "Arid3a(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = TRUE, repel = TRUE
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = plot_path, width = 4, height = 4, units = "in", dpi = 150)


# auc_mtx is exported from pySCENIC
# this aut_mtx is generated from the analysis in Epi-PC and Epi-PTC alone
auc_mtx <- read_tsv(file = file.path("int_data", "pyscenic_nig-sc", "epi_PC_PTC", "auc_mtx.txt")) %>% 
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
DimPlot(seu, group.by = "celltype_2", label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype_2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, "regulon_celltype.xlsx"))

Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
markers <- FindAllMarkers(seu, only.pos = TRUE)

DefaultAssay(seu) <- "scenicAUC"
add_feat <- "Foxm1(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = TRUE, repel = TRUE
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = plot_path, width = 4, height = 4, units = "in", dpi = 150)
