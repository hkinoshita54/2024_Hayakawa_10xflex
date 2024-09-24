analysis_step <- "080_combine_annotation"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper functions ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}
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

# Add annotation ----
## Load annotated seurat objects
epi <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
imm <- readRDS(file.path("RDSfiles", "seu_030.1_imm.RDS"))
str <- readRDS(file.path("RDSfiles", "seu_040.1_str.RDS"))
epi_PC_PTC <- readRDS(file.path("RDSfiles", "seu_021.1_epi_PC&PTC.RDS"))
endo <- readRDS(file.path("RDSfiles", "seu_041.1_endothelial.RDS"))
fibro <- readRDS(file.path("RDSfiles", "seu_042.1_fibroblast.RDS"))
seu <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))

## Add annotation as character
seu$cellgroup <- "Epi."
seu$cellgroup[colnames(imm)] <- "Imm."
seu$cellgroup[colnames(str)] <- "Str."

seu$celltype <- ""
seu$celltype[colnames(epi)] <- as.character(epi$celltype)
seu$celltype[colnames(imm)] <- as.character(imm$celltype)
seu$celltype[colnames(str)] <- as.character(str$celltype)

seu$celltype_2 <- seu$celltype
seu$celltype_2[colnames(epi_PC_PTC)] <- as.character(epi_PC_PTC$celltype_2)
seu$celltype_2[colnames(endo)] <- as.character(endo$celltype_2)
seu$celltype_2[colnames(fibro)] <- as.character(fibro$celltype_2)

seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Imm.", "Str."))
seu$celltype <- factor(seu$celltype, levels = c(levels(epi$celltype), levels(imm$celltype), levels(str$celltype)))
seu$celltype_2 <- factor(seu$celltype_2, 
                         levels = c(levels(epi$celltype)[-10],
                                    levels(epi_PC_PTC$celltype_2),
                                    levels(imm$celltype),
                                    levels(fibro$celltype_2),
                                    levels(endo$celltype_2),
                                    levels(str$celltype)[c(2, 5:10)]
                                    )
                         )
seu <- seu[, !is.na(seu$celltype_2)]

# Clustering ----
seu <- cluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150) 
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "cellgroup") & NoAxes()
ggsave("cellgroup.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype", cols = "polychrome") + NoAxes()  +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("celltype.png", path = plot_path, width = 6, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype_2", cols = "polychrome") + NoAxes()  +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("celltype_2.png", path = plot_path, width = 7, height = 3, units = "in", dpi = 150)

# Dot plot
features <- c("Epcam", "Cldn18", "Ptprc", "Arhgap45", "Dcn", "Sparc")
DotPlot(seu, group.by = "cellgroup", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 4, height = 3.7, units = "in", dpi = 150)

# Save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_080_combined.RDS"))
