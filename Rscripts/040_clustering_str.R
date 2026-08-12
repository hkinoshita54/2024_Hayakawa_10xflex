# Continued from 010.1_clustering.R
# Without integration
analysis_step <- "040_clustering_str"

# load packages ----
library(tidyverse)
library(Seurat)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper function ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
pal <- readRDS("RDSfiles/color_palette.RDS")
load(file = "RDSfiles/cellcycle_genes_mouse.Rdata")
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))
seu <- seu_all[, str_names]
seu$cellgroup <- "Str."

# Clustering ----
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- RunPCA(seu, npcs = 50)
seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = 1)
seu <- RunUMAP(seu, dims = 1:50)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

# adjust resolution
# seu <- FindClusters(seu, resolution = 3, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) & NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)


# Check markers by feature plots
features <- readLines("aux_data/gene_set/annotation/03_str_markers.txt")
sapply(features, save_fp, seu, fp_path)

# add_feat <- "Pgc"
# FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
# ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype1 annotation
seu$celltype1 <-""
seu$celltype1[seu$seurat_clusters %in% c(4,7,18)] <- "Fibro."
seu$celltype1[seu$seurat_clusters %in% c(0,1,5,6,14,15,17,27)] <- "Fibro."
seu$celltype1[seu$seurat_clusters %in% c(2,8,9,19,21)] <- "EC"
seu$celltype1[seu$seurat_clusters %in% c(3,11,12,13,20,24)] <- "Myo."
seu$celltype1[seu$seurat_clusters %in% c(10)] <- "Peri."
seu$celltype1[seu$seurat_clusters %in% c(22)] <- "Adipo."
seu$celltype1[seu$seurat_clusters %in% c(16)] <- "Glial"
seu$celltype1[seu$seurat_clusters %in% c(25,26)] <- "Neural"
seu$celltype1[seu$seurat_clusters %in% c(23)] <- "Meso."
seu$celltype1 <- factor(seu$celltype1, levels = c("Fibro.", "EC", "Myo.", "Peri.", "Adipo.", "Glial", "Neural", "Meso."))
DimPlot(seu, group.by = "celltype1", label = TRUE, repel = TRUE, cols = pal$stromal) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype1.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)


# save the cell barcodes of each celltype1
Fib_names <- colnames(seu)[seu$celltype1 ==  "Fibro."]
EC_names <- colnames(seu)[seu$celltype1 == "EC"]
Myo_names <- colnames(seu)[seu$celltype1 == "Myo."]
Peri_names <- colnames(seu)[seu$celltype1 == "Peri."]
Adipo_names <- colnames(seu)[seu$celltype1 == "Adipo."]
Glial_names <- colnames(seu)[seu$celltype1 == "Glial"]
Neural_names <- colnames(seu)[seu$celltype1 == "Neural"]
Meso_names <- colnames(seu)[seu$celltype1 == "Meso."]
save(Fib_names, EC_names, Myo_names, Peri_names, Adipo_names, Glial_names, Neural_names, Meso_names,
     file = file.path("RDSfiles", "str_celltype1_names.Rdata"))

# Dot plot
Idents(seu) <- "celltype1"
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))
features <- c("Pdgfra","Col1a1",
              "Pecam1","Egfl7",
              "Myh11","Acta2",
              "Rgs5","Notch3",
              "Adipoq","Plin1",
              "Sox10","Plp1",
              "Kit","Etv1",
              "Wt1","Msln")
DotPlot(seu, group.by = "celltype1", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6.5, height = 4.5, units = "in", dpi = 150)
