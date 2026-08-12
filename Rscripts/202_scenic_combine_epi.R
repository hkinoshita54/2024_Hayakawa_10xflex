# Continued from 020_clustering_epi.R and pySCENIC output
analysis_step <- "202_scenic_combine_epi"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(ggrepel)
source("Rscripts/theme_panel.R")

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path("plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# load data ----
seu <- readRDS("RDSfiles/seu_053_epi_combined.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

# auc_mtx is exported from pySCENIC
auc_mtx <- read_tsv(file = "pyscenic_docker/053_combine_epi/hvg8000_alone/adj_top200/auc_mtx.txt") %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[,Cells(seu)]
seu[["scenicAUC"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "scenicAUC"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu_cp <- seu
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype2", cols = pal$epi, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend()
ggsave("celltype2.pdf", path = plot_path, width = 20, height = 30, units = "mm")

DimPlot(seu, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("regulon_celltype.xlsx")))

# check RNA markers of scenic identified clusters
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("scenic_cluster_on_Seurat_UMAP.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("markers_scenic_clusters.xlsx")))

# feature plot
add_feat <- "Rel(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q50", max.cutoff = "q75"
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype2"
# saveRDS(seu, file = file.path("RDSfiles", "seu_202_scenic_combine_epi.RDS"))

# edit rss.csv
rss <- read_csv("pyscenic_docker/053_combine_epi/rss.csv")
rss[,1] <- sapply(rss[,1], gsub, pattern = "_", replacement = "")
colnames(rss)[1] <- "Regulon"
openxlsx2::write_xlsx(rss, file.path(res_path, paste0("rss.xlsx")))
