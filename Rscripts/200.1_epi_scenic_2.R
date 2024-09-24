# Continued from 020.1_clustering_epi_2.R and pySCENIC output (nig-sc)
analysis_step <- "200.1_epi_scenic_2"

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
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))

# auc_mtx is exported from pySCENIC (5th iteration)
auc_mtx <- read_tsv(file = file.path("int_data", "pyscenic_nig-sc", "epi_2", "auc_mtx.txt")) %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
seu[["scenicAUC"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "scenicAUC"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)

DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype"
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
fp_path <- file.path("plot", analysis_step, "feature_plot")
fs::dir_create(fp_path)
DefaultAssay(seu) <-"scenicAUC"

add_feat <- "Sox4(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75"
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype"
saveRDS(seu, file = file.path("RDSfiles", "seu_200.1_epi_scenic.RDS"))

# edit rss.csv
rss <- read_csv("int_data/pyscenic_nig-sc/epi_2/rss.csv")
rss[,1] <- sapply(rss[,1], gsub, pattern = "_", replacement = "")
openxlsx2::write_xlsx(rss, file.path(res_path, paste0("rss.xlsx")))
