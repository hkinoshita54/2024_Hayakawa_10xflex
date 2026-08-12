# Continued from 020_clustering_epi.R
# Without integration
analysis_step <- "021_clustering_epi_PC&PTC"

# load packages ----
library(tidyverse)
library(Seurat)
library(speckle)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
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

# Clustering ----

## 1. Initial clustering
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))
seu <- subset(seu, subset = celltype1 %in% c("Epi-PC&PTC") & orig.ident %in% c("PC", "PTC"))
seu <- recluster(seu, npcs = 50, res = 0.5)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) + NoAxes()
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.5.csv"))

## 2. Adjust resolution
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res1.csv"))

## 3. remove contamination
seu <- subset(seu, subset = seurat_clusters %in% c(10,11,12), invert = TRUE)
seu <- recluster(seu, npcs = 50, res = 0.5)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) + NoAxes()
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.5_2.csv"))

## 4. Adjust resolution again
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res1_2.csv"))

## 5. Decided to use this clustering and save plots
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4) & NoAxes()
ggsave("sample.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

# Add celltyp3 annotation and save the Seurat object ----
seu$celltype3 <- ""
seu$celltype3[seu$seurat_clusters %in% c(6)] <- "Dys.1"
seu$celltype3[seu$seurat_clusters %in% c(1,7)] <- "Dys.2"
seu$celltype3[seu$seurat_clusters %in% c(4)] <- "Dys.3"
seu$celltype3[seu$seurat_clusters %in% c(2,9)] <- "Trans."
seu$celltype3[seu$seurat_clusters %in% c(8)] <- "Ca.1"
seu$celltype3[seu$seurat_clusters %in% c(5)] <- "Ca.2"
seu$celltype3[seu$seurat_clusters %in% c(0)] <- "Ca.3"
seu$celltype3[seu$seurat_clusters %in% c(3)] <- "Ca.4"
seu$celltype3 <- factor(seu$celltype3, levels = c("Dys.1", "Dys.2", "Dys.3", "Trans.", 
                                                  "Ca.1", "Ca.2", "Ca.3", "Ca.4"))
DimPlot(seu, group.by = "celltype3", cols = pal$tum, label = TRUE, repel = TRUE) + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype3.png", path = plot_path, width = 3.5, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype3"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

seu$RNA_snn_res.2 <- NULL
saveRDS(seu, file = file.path("RDSfiles", "seu_021_tumor.RDS"))

features <- c(
  "Ankle1", "Top2a", "Cit", "Kif18b", "Ect2",
  "Onecut2", "Tmprss2", "Muc3a", "Ern2", "Kcne3",
  "Gabarap", "Eef1a1", "Wfdc18", "Gstm1", "Actr2",
  "Mki67", "Ccna2", "Aurkb", "Plk1", "Ube2c",
  "Lgr5", "Lrg1", "Cd177", "Reg3g", "Gkn3",
  "Cldn4", "Krt8", "Tnfrsf12a", "Lamc2", "Psca"
)

seu_dot <- subset(seu, subset = celltype3 %in% c("Trans.", "Ca.4"), invert = TRUE)
DotPlot(seu_dot, group.by = "celltype3", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 9, height = 4, units = "in", dpi = 150)

add_feat <- "Top2a"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Cell numbers and proportion plot ----
seu <- readRDS("RDSfiles/seu_021_tumor.RDS")
fraction <- "tumor"

nclust <- length(levels(seu$celltype3))
props <- getTransformedProps(clusters = seu$celltype3, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$tum) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 2.3, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$tum) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 2.3, height = 3.2, units = "in", dpi = 150)
