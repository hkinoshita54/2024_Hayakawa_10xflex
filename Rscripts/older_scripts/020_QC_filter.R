####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")


####################
# load data and filter by doublet calls, nFeature_RNA and percent.mt ----
path = "plots/020_QC_filter"
dir.create(path)
seu <- readRDS(file = "RDSfiles/seu_01_unfiltered_DoubletFinder.RDS")
seu@meta.data <- seu@meta.data[,c(1:4,8)]
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_unfiltered.png", path = path, width = 9, height = 3, units = "in", dpi = 150)
seu <- subset(seu, subset = doublet_finder == "Singlet")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_singlet.png", path = path, width = 9, height = 3, units = "in", dpi = 150)

# filter arbitrarily
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
# seu <- subset(seu, subset = Ptprc > 0, invert = TRUE)
# seu <- subset(seu, subset = Pecam1 > 0, invert = TRUE)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_RNA_mt_filtered.png", path = path, width = 9, height = 3, units = "in", dpi = 150)


saveRDS(seu, file = "RDSfiles/seu_02_filtered.RDS")
