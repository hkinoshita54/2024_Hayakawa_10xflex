# doing this AFTER 053
analysis_step <- "052_combine_all"

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
source("Rscripts/theme_panel.R")
cluster = function(seu, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
} 

recluster = function(seu, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
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

pal <- readRDS("RDSfiles/color_palette.RDS")

# Add annotation ----
## Load seurat objects, all cells and annotated cells
seu_all <- readRDS("RDSfiles/seu_010_filt.RDS")

epi <- readRDS("RDSfiles/seu_053_epi_combined.RDS")
epi_names <- Cells(epi)

imm <- readRDS("RDSfiles/seu_050_imm_combined.RDS")
imm_names <- Cells(imm)

str <- readRDS("RDSfiles/seu_051_str_combined.RDS")
str_names <- Cells(str)

seu <- subset(seu_all, cells = c(epi_names, imm_names, str_names))

## Add annotation
### cellgroup
seu$cellgroup <- ""
seu$cellgroup[epi_names] <- "Epi."
seu$cellgroup[imm_names] <- "Imm."
seu$cellgroup[str_names] <- "Str."
seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Imm.", "Str."))

### celltype1
seu$celltype1 <- ""
seu$celltype1[epi_names] <- as.character(epi$celltype1)
seu$celltype1[imm_names] <- as.character(imm$celltype1)
seu$celltype1[str_names] <- as.character(str$celltype1)
epi_levels <- levels(epi$celltype1)
imm_levels <- levels(imm$celltype1)
str_levels <- levels(str$celltype1)
seu$celltype1 <- factor(seu$celltype1, levels = c(epi_levels, imm_levels, str_levels))

### celltype2
seu$celltype2 <- ""
seu$celltype2[epi_names] <- as.character(epi$celltype2)
seu$celltype2[imm_names] <- as.character(imm$celltype2)
seu$celltype2[str_names] <- as.character(str$celltype2)
epi_levels <- levels(epi$celltype2)
imm_levels <- levels(imm$celltype2)
str_levels <- levels(str$celltype2)
seu$celltype2 <- factor(seu$celltype2, levels = c(epi_levels, imm_levels, str_levels))

# save minimal seurat object ----
saveRDS(seu, "RDSfiles/seu_052_all_combined.RDS")


# Clustering & plots ----
seu <- cluster(seu, npcs = 50, res = 1)
saveRDS(seu, "RDSfiles/seu_052_all_combined.RDS")

seu <- readRDS("RDSfiles/seu_052_all_combined.RDS")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "cellgroup", cols = pal$cellgroup, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes()
ggsave("cellgroup.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype1", label = TRUE, repel = TRUE, 
        cols = c(pal$epi[1:11], unname(pal$immune_group), pal$stromal)) &
  NoAxes() & NoLegend()
ggsave("celltype1.png", path = plot_path, width = 3, height = 3, units = "in", dpi = 150)

DimPlot(seu, group.by = "celltype2", label = TRUE, repel = TRUE, 
        cols = c(pal$epi, 
                 pal$T_cells, pal$B_cells, pal$myeloid,
                 pal$fibro_sub, unname(pal$endothelial), pal$stromal[3:8])) & 
  NoAxes() & 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("celltype2.png", path = plot_path, width = 10, height = 5, units = "in", dpi = 150)

Idents(seu) <- "cellgroup"
features <- c("Epcam", "Cldn18", "Ptprc", "Arhgap45", "Dcn", "Sparc")
features <- c("Cldn18", "Arhgap45","Sparc")
DotPlot(seu, group.by = "cellgroup", features = features, dot.scale = 2.3) + 
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
ggsave("dotplot2.pdf", path = plot_path, width = 28, height = 35, units = "mm")

# Cell numbers and proportion plot ----
fraction <- "all"

nclust <- length(levels(seu$cellgroup))
props <- getTransformedProps(clusters = seu$cellgroup, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$cellgroup) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 35, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$cellgroup) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 35, units = "mm")

# Feature plot for EDA ----
add_feat <- "Tgfb1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)
