# Continued from 010_QC.R
# Without integration
analysis_step <- "020_clustering_epi"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))

# Clustering
seu <- seu_all[, epi_names]
seu$cellgroup <- "Epi."
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "annotation", "01_epi_markers.txt"))

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
         width = 4, height = 4, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
# seu <- JoinLayers(seu)
Idents(seu) <- "seurat_clusters"
markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Tgfbr2"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = TRUE, repel = TRUE
            ) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)
rm(add_feat)

# Add celltype annotation and save the Seurat object
seu <- subset(seu, idents = 21, invert = TRUE)
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(14)] <- "Prog."
seu$celltype[seu$seurat_clusters %in% c(6,8)] <- "Pit"
seu$celltype[seu$seurat_clusters %in% c(1,9)] <- "Neck"
seu$celltype[seu$seurat_clusters %in% c(2,4)] <- "Chief"
seu$celltype[seu$seurat_clusters %in% c(5,16)] <- "Pariet."
seu$celltype[seu$seurat_clusters %in% c(17,19)] <- "EEC"
seu$celltype[seu$seurat_clusters %in% c(20)] <- "Tuft"
seu$celltype[seu$seurat_clusters %in% c(15,18)] <- "Squam."
seu$celltype[seu$seurat_clusters %in% c(11)] <- "Epi-PT"
seu$celltype[seu$seurat_clusters %in% c(0,13)] <- "Epi-PC"
seu$celltype[seu$seurat_clusters %in% c(3,7,10,12)] <- "Epi-PTC"
seu$celltype <- factor(seu$celltype, levels = c("Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Tuft",
                                                "Squam.", "Epi-PT", "Epi-PC", "Epi-PTC"))
DimPlot(seu, group.by = "celltype", cols = "alphabet", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype"
markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("Mki67", "Muc5ac", "Muc6", "Cblif", "Atp4b", "Chga", "Dclk1", "Krt5",
              "Ces2a", "Cym", "Lrg1")
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))
