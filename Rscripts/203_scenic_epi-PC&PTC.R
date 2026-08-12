# Continued from 021_clustering_epi-PC&PTC.R and pySCENIC output
analysis_step <- "203_scenic_epi-PC&PTC"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(ggrepel)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path("plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# load data ----
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

# auc_mtx is exported from pySCENIC (5th iteration)
auc_mtx <- read_tsv(file = "pyscenic_docker/023_tumor_ssGSEA/auc_mtx.txt") %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[,Cells(seu)]
seu[["reg"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "reg"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu_cp <- seu
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 4) + NoAxes()
ggsave("cluster.pdf", path = plot_path, width = 3.3, height = 3, units = "in")
DimPlot(seu, group.by = "celltype3", cols = pal$tum, raster = TRUE, pt.size = 4) + NoAxes()
ggsave("celltype3.pdf", path = plot_path, width = 3.5, height = 3, units = "in")
DimPlot(seu, group.by = "orig.ident", raster = TRUE, pt.size = 4) + NoAxes()
ggsave("sample.pdf", path = plot_path, width = 3.5, height = 3, units = "in")

Idents(seu) <- "celltype3"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("regulon_celltype.xlsx")))

# check RNA markers of scenic identified clusters
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 4) + NoAxes()
ggsave("scenic_cluster_on_Seurat_UMAP.pdf", path = plot_path, width = 3.3, height = 3, units = "in")

markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("markers_scenic_clusters.xlsx")))

# feature plot and vln ----
add_feat <- "Tead4(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q50", max.cutoff = "q75"
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

VlnPlot(seu, group.by = "celltype3", features = add_feat, cols = pal$tum, pt.size = 0) & NoLegend()
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 3, height = 3, units = "in")

# save seurat object ----
Idents(seu) <- "celltype3"
saveRDS(seu, file = file.path("RDSfiles", "seu_203_scenic_epi-PC&PTC.RDS"))

# edit rss.csv ---
rss <- read_csv("pyscenic_docker/023_tumor_ssGSEA/rss.csv")
rss[,1] <- sapply(rss[,1], gsub, pattern = "_", replacement = "")
colnames(rss)[1] <- "Regulon"
openxlsx2::write_xlsx(rss, file.path(res_path, paste0("rss.xlsx")))
