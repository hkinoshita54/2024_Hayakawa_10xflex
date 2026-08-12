# Continued from 010.1_clustering.R
# Without integration
analysis_step <- "020_clustering_epi"

# load packages ----
library(tidyverse)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper function ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

# recluster_sct = function(seu_obj, npcs, res, n.neighbors, min.dist){
#   seu[["RNA"]]$scale.data <- NULL
#   seu[["RNA"]]$data <- NULL
#   seu <- SCTransform(seu, vst.flavor = "v2", verbose = F, variable.features.n = 4000)
#   seu <- RunPCA(seu, npcs = npcs)
#   seu <- FindNeighbors(seu, dims = 1:npcs)
#   seu <- FindClusters(seu, resolution = res)
#   seu <- RunUMAP(seu, dims = 1:npcs, n.neighbors = n.neighbors, min.dist = min.dist)
#   return(seu)
# }

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Load files ----
pal <- readRDS("RDSfiles/color_palette.RDS")
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))
seu <- seu_all[, epi_names]
seu$cellgroup <- "Epi."

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE)
seu_cp <- seu

seu <- subset(seu, subset = seurat_clusters %in% 26:29, invert = TRUE) # pancreas, neuron, enterocyte, etc.
seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
markers <- FindAllMarkers(seu, only.pos = TRUE)
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) & NoAxes()
ggsave("sample.png", path = plot_path, width = 3.3, height = 3, units = "in", dpi = 150)

# # see how they are clustered with SCTransform
# seu <- recluster_sct(seu, npcs = 20, res = 2, n.neighbors = 50, min.dist = 0.4)
# DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
# markers <- FindAllMarkers(seu, only.pos = TRUE)
# ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
# DimPlot(seu, group.by = "orig.ident") + NoAxes()
# ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots ----
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "01_epi_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Lrg1", "Cd38", "Cldn4", "Msln", "Cxcl5", "Pigr")
sapply(add_feat, save_fp, seu, fp_path)
seu_cp2 <- seu

# seu <- recluster(seu, npcs = 50, res = 3)
# DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
# ggsave("cluster_res2.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

# Add celltype1 annotation ----
seu$celltype1 <- ""
seu$celltype1[seu$seurat_clusters %in% c(9)] <- "Stem"
seu$celltype1[seu$seurat_clusters %in% c(13,17,33)] <- "Prog."
seu$celltype1[seu$seurat_clusters %in% c(10,12,24)] <- "Pit"
seu$celltype1[seu$seurat_clusters %in% c(2,4)] <- "Neck"
seu$celltype1[seu$seurat_clusters %in% c(0,15,16,22)] <- "Chief"
seu$celltype1[seu$seurat_clusters %in% c(5,18,25)] <- "Pariet."
seu$celltype1[seu$seurat_clusters %in% c(26,29,31)] <- "EEC"
seu$celltype1[seu$seurat_clusters %in% c(34)] <- "Tuft"
seu$celltype1[seu$seurat_clusters %in% c(19,20)] <- "Squam."
seu$celltype1[seu$seurat_clusters %in% c(21,23)] <- "Pit-PT"
seu$celltype1[seu$seurat_clusters %in% c(1,3,6:8,11,14,27,28,30,32)] <- "Epi-PC&PTC"
seu$celltype1 <- factor(seu$celltype1, levels = c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                  "Squam.", "Pit-PT", "Epi-PC&PTC"))
DimPlot(seu, group.by = "celltype1", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype1.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Add celltype2 annotation ----
tum_names <- colnames(seu)[seu$celltype1 %in% c("Epi-PC&PTC")]
seu$celltype2 <- as.character(seu$celltype1)
seu$celltype2[tum_names] <- paste0("Epi-", seu$orig.ident[tum_names])
seu <- subset(seu, subset = celltype2 %in% c("Epi-WT", "Epi-PT"), invert = TRUE)
seu$celltype2 <- factor(seu$celltype2, levels = c("Stem", "Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                  "Squam.", "Pit-PT", "Epi-PC", "Epi-PTC"))
DimPlot(seu, group.by = "celltype2", cols = pal$epi, label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot ----
Idents(seu) <- "celltype2"

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

features <- c("Mki67", "Muc5ac", "Muc6", "Cblif", "Atp4b", "Chga", "Dclk1", "Krt5", "Ces2a", "Pigr")
DotPlot(seu, group.by = "celltype1", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

features <- c(
  "Mki67","Top2a",
  "Tff2","Muc5ac",
  "Gkn1","Mucl3",
  "Muc6","Pga5",
  "Pgc","Clps",
  "Atp4a","Atp4b",
  "Chga","Chgb",
  "Pou2f3","Trpm5",
  "Krt5","Trp63",
  "Sctr","Ces2a",
  "Pigr","Muc4",
  "Lrg1","Cxcl5"
)

DotPlot(seu, group.by = "celltype2", features = features) + RotatedAxis()
ggsave("dotplot_2.png", path = plot_path, width = 8, height = 4, units = "in", dpi = 150)

# save RDS file ----
seu$RNA_snn_res.1 <- NULL
saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))

# Cell numbers and proportion plot ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")
fraction <- "epi"

nclust <- length(levels(seu$celltype2))
props <- getTransformedProps(clusters = seu$celltype2, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$epi) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 2.7, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$epi) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 2.7, height = 3.2, units = "in", dpi = 150)
