####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")


####################
# set conditions ----
npcs = 15
path = paste0("plots/410_T1PTC_exp3T1P_lognorm_harmony_", "npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_410_T1PTC_exp3T1P_lognorm_harmony_npc", as.character(npcs), ".RDS")


####################
# load data and filter by doublet calls, nFeature_RNA and percent.mt ----
seu <- readRDS(file = "RDSfiles/seu_02_filtered.RDS")
seu <- subset(seu, subset = orig.ident == "T1PTC")
seu2 <- readRDS(file = "RDSfiles/prev_exp/seu_T1P_22_filtered.RDS")
seu <- merge(x = seu, y = seu2)


####################
# lognormalization with harmony integration ----
seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay5")
seu[["RNA"]]$data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident) 
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seucp <- seu
# seu <- seucp
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, reduction = "harmony", resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
# DimPlot(seu, group.by = "batch") + NoAxes()
# ggsave("batch.png", path = path, width = 5, height = 5, units = "in", dpi = 150)


####################
#### feature plots ####
files <- list.files(path = "gene_set/annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


####################
#### violin plots ----
# for(i in 1:length(features)){
#   tryCatch({
#     p <- VlnPlot(seu, features = features[i])
#     ggsave(paste0("vln_", as.character(features[i]), ".png"), plot = p,
#            path = path,
#            width = 9, height = 3, units = "in", dpi = 150)
#   }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
# }
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0)
# ggsave("vln_RNA_mt.png", path = path, width = 12, height = 9, units = "in", dpi = 150)


saveRDS(seu, file = RDSfile)
