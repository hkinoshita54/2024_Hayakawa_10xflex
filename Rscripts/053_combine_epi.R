analysis_step <- "053_combine_epi"

# load packages ----
library(tidyverse)
library(Seurat)
library(speckle)
source("Rscripts/theme_panel.R")

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
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 4) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = fp_path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

pal <- readRDS("RDSfiles/color_palette.RDS")

# Add annotation ----
## Load seurat objects, all cells and annotated immune cells

# all the epithelial cells from 020
epi_all <- readRDS("RDSfiles/seu_020_epi.RDS")
tum <- subset(epi_all, subset = celltype1 %in% c("Epi-PC&PTC") & orig.ident %in% c("PC", "PTC")) #tumor subset from epi_all
ncol(tum)

# filtered tum_filt from 021
tum_filt <- readRDS("RDSfiles/seu_021_tumor.RDS") 
ncol(tum_filt)

# cells to be removed from 
filt_names <- setdiff(Cells(tum), Cells(tum_filt))
length(filt_names)

# remove filtered cells in 021
epi <- subset(epi_all, cells = filt_names, invert = TRUE)
ncol(epi)
names(epi[[]])

## Add annotation
### tumor subsets
epi$celltype3 <- as.character(epi$celltype1)
epi$celltype3[Cells(tum_filt)] <- as.character(tum_filt$celltype3)
epi$celltype3 <- factor(epi$celltype3, 
                        levels = c(levels(epi$celltype1)[-11], levels(tum_filt$celltype3)))

# plots ----
seu <- epi
# seu <- recluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype1", cols = pal$epi, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype1.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype2", cols = pal$epi, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype2.pdf", path = plot_path, width = 45, height = 30, units = "mm")

Idents(seu) <- "celltype1"
features <- c(
  "Top2a",
  "Tff2",
  "Muc5ac",
  "Muc6",
  "Cblif",
  "Atp4a",
  "Chga",
  "Trpm5",
  "Krt5",
  "Ces2a",
  "Pigr"
)
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  )+
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot2.pdf", path = plot_path, width = 50, height = 40, units = "mm")

# save seurat object ----
# seu <- DietSeurat(seu)
# seu[["RNA"]]$scale.data <- NULL
# seu[["RNA"]]$data <- NULL
saveRDS(seu, "RDSfiles/seu_053_epi_combined.RDS")

# Cell numbers and proportion plot ----
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
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 42, height = 36, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$epi) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 28, height = 36, units = "mm")

# EDA ----
seu <- readRDS("RDSfiles/seu_053_epi_combined.RDS")

add_feat <- c("Spp1")
sapply(add_feat, save_fp, seu, fp_path)

VlnPlot(seu, group.by = "celltype2", features = add_feat, cols = pal$epi, pt.size = 0) &
  theme_panel() & NoLegend() &
  labs(x = NULL, y = "Expression") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 45, height = 30, units = "mm")

