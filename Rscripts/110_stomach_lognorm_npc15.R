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
path = paste0("plots/110_stomach_lognorm_", "npc", as.character(npcs))
dir.create(path)
RDSfile = paste0("RDSfiles/seu_110_stomach_lognorm_npc", as.character(npcs), ".RDS")


####################
# load data and filter by doublet calls, nFeature_RNA and percent.mt ----
seu <- readRDS(file = "RDSfiles/seu_02_filtered.RDS")
seu <- subset(seu, subset = exp_group == "stomach")


####################
# cluster without integration ----
seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay5")
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
# DimPlot(seu, group.by = "exp_group") + NoAxes()
# ggsave("group.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
# DimPlot(seu, split.by = "exp_group") + NoAxes()
# ggsave("group_split.png", path = path, width = 8, height = 4, units = "in", dpi = 150)


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


####################
#### additional feature plots ####
seu <- readRDS(RDSfile)
files <- list.files(path = "gene_set/additional/", full.names = TRUE)
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
