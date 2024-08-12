####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")

####################
# set conditions ----
description = "stomach_lognorm_harmony"
file_no = "120"
npcs = 15
path = file.path("plots", paste0(file_no, "_", description, "_npc", as.character(npcs)) )
dir.create(path = path, recursive = TRUE)
RDSfile = paste0("RDSfiles/", "seu_", file_no, "_", description, "_npc", as.character(npcs), ".RDS")


####################
# load data and merge, then split by id ----
seu <- readRDS(file = "RDSfiles/seu_02_filtered.RDS")
seu <- subset(seu, subset = exp_group == "stomach")


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
for(i in 1:length(features)){
  tryCatch({
    p <- VlnPlot(seu, features = features[i])
    ggsave(paste0("vln_", as.character(features[i]), ".png"), plot = p,
           path = path,
           width = 9, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0)
ggsave("vln_RNA_mt.png", path = path, width = 12, height = 9, units = "in", dpi = 150)


saveRDS(seu, file = RDSfile)
