# Continued from 041.1_clustering_endothelial_2.R
# combine TumEC1&2 together to do the downstream analyses
analysis_step <- "041.2_clustering_endothelial_3"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(escape)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper function ----
recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
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

# Load data ----
seu <- readRDS(file.path("RDSfiles", "seu_041.1_endothelial.RDS"))
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()

# rename clusters
seu$celltype_3 <- seu$celltype_2   # stash original celltype_2 in celltype_3
seu$celltype_2 <- as.character(seu$celltype_2)
seu$celltype_2[seu$celltype_2 %in% c("TumEC1", "TumEC2")] <- "TumEC"
seu$celltype_2[seu$celltype_2 %in% c("TumEC3")] <- "Prolif.EC"
seu$celltype_2 <- factor(seu$celltype_2, levels = c("ArtEC", "VenEC", "CapEC", "TumEC", "Prolif.EC", "LEC", "TumLEC"))
DimPlot(seu, group.by = "celltype_2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
Idents(seu) <- "celltype_2"

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- c("Tnfrsf9")
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Il4", "Cd28", "Cadm3", "Spon2", "Icosl", "Lta", "Btla", "Tgm2", "Ccl2")
add_feat <- c("Lrg1", "Eng", "Cd276", "Tnfrsf9", "Sele", "Selp", "Cxcl12", "Cd274", "Gpnmb", "Pdcd1")
sapply(add_feat, save_fp, seu, fp_path)


# Check and save markers
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1000) %>%
  ungroup() -> markers
openxlsx2::write_xlsx(markers, file.path(res_path, "markers.xlsx"))

# dotplot
features = c("Sema3g", "Ackr1", "Rgcc", "Kit", "Top2a", "Prox1")
DotPlot(seu, group.by = "celltype_3", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5, height = 4, units = "in", dpi = 150)

# save data
saveRDS(seu, file = file.path("RDSfiles", "seu_041.1_endothelial.RDS"))


# ss enrichment by escape ----
H <- getGeneSets(species = "Mus musculus", library = "H")
names(H) <- gsub("HALLMARK-", "", names(H))
names(H) <- gsub("-", "_", names(H))

# ssGSEA ----
seu_all <- seu
seu <- subset(seu_all, subset = celltype_3 %in% c("ArtEC", "VenEC", "CapEC", "TumEC", "Prolif.EC"))
seu <- runEscape(seu, method = "ssGSEA", gene.sets = H,
                 groups = 5000, min.size = 15, new.assay.name = "ssGSEA_H",
                 BPPARAM = SnowParam(workers = 2))
seu <- performNormalization(seu, assay = "ssGSEA_H", gene.sets = H, 
                            scale.factor = seu$nFeature_RNA)

gs_H <- FindAllMarkers(seu, assay = "ssGSEA_H_normalized", min.pct = 0, logfc.threshold = 0)
openxlsx2::write_xlsx(gs_H, file.path(res_path, "ssGSEA_H.xlsx"))

# plot heatmap
heatmapEnrichment(seu, group.by = "ident", assay = "ssGSEA_H_normalized", gene.set.use = "all",
                  scale = TRUE, cluster.rows = TRUE, cluster.columns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "right", legend.direction = "vertical")
ggsave("GSEA_H_hm.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)


# ssGSEA without Prolif.EC ----
seu <- subset(seu_all, subset = celltype_3 %in% c("ArtEC", "VenEC", "CapEC", "TumEC"))
seu <- runEscape(seu, method = "ssGSEA", gene.sets = H,
                 groups = 5000, min.size = 15, new.assay.name = "ssGSEA_H",
                 BPPARAM = SnowParam(workers = 2))
seu <- performNormalization(seu, assay = "ssGSEA_H", gene.sets = H, 
                            scale.factor = seu$nFeature_RNA)

gs_H <- FindAllMarkers(seu, assay = "ssGSEA_H_normalized", min.pct = 0, logfc.threshold = 0)
openxlsx2::write_xlsx(gs_H, file.path(res_path, "NoProlifEC_ssGSEA_H.xlsx"))

# plot heatmap
heatmapEnrichment(seu, group.by = "ident", assay = "ssGSEA_H_normalized", gene.set.use = "all",
                  scale = TRUE, cluster.rows = TRUE, cluster.columns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "right", legend.direction = "vertical")
ggsave("NoProlifEC_GSEA_H_hm.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)