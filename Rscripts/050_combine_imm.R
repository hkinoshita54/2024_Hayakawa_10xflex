analysis_step <- "050_combine_imm"

# load packages ----
library(tidyverse)
library(Seurat)
library(speckle)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Helper functions ----
source("Rscripts/theme_panel.R")
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

pal <- readRDS("RDSfiles/color_palette.RDS")

# Add annotation ----
## Load seurat objects, all cells and annotated immune cells
seu_all <- readRDS("RDSfiles/seu_010_filt.RDS")
tcell <- readRDS("RDSfiles/seu_031_tcell.RDS")
bcell <- readRDS("RDSfiles/seu_032_bcell.RDS")
myeloid <- readRDS("RDSfiles/seu_033_myeloid.RDS")

## subset to all the three immune cells
t_names <- Cells(tcell)
b_names <- Cells(bcell)
m_names <- Cells(myeloid)
imm_names <- c(t_names, b_names, m_names)
seu <- subset(seu_all, cells = imm_names)

## Add annotation
### immune groups
seu$celltype1 <- "Mye."
seu$celltype1[b_names] <- "Bcell"
seu$celltype1[t_names] <- "Tcell"
seu$celltype1 <- factor(seu$celltype1, levels = c("Tcell", "Bcell", "Mye."))

### annotation within tcell, bcell, or myeloid
seu$celltype2 <- ""
seu$celltype2[t_names] <- as.character(tcell$celltype2)
seu$celltype2[b_names] <- as.character(bcell$celltype2)
seu$celltype2[m_names] <- as.character(myeloid$celltype2)
t_levels <- levels(tcell$celltype2)
b_levels <- levels(bcell$celltype2)
m_levels <- levels(myeloid$celltype2)
all_levels <- c(t_levels, b_levels, m_levels)
seu$celltype2 <- factor(seu$celltype2, levels = all_levels)

# save minimal seurat object ----
saveRDS(seu, "RDSfiles/seu_050_imm_combined.RDS")

# Clustering & plots ----
seu <- cluster(seu, npcs = 50, res = 1)
saveRDS(seu, "RDSfiles/seu_050_imm_combined.RDS")

seu <- readRDS("RDSfiles/seu_050_imm_combined.RDS")

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype1", cols = pal$immune_group, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype1")
ggsave("celltype1.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype2", 
        cols = c(pal$T_cells, pal$B_cells, pal$myeloid), 
        raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype2")
ggsave("celltype2.pdf", path = plot_path, , width = 60, height = 45, units = "mm")

Idents(seu) <- "celltype1"
features <- c("Cd3e", "Cd3g", "Cd19", "Pax5", "Csf2ra", "Tyrobp")
DotPlot(seu, group.by = "celltype1", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 45, height = 45, units = "mm")

# Cell numbers and proportion plot ----
fraction <- "imm"

nclust <- length(levels(seu$celltype1))
props <- getTransformedProps(clusters = seu$celltype1, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$immune_group) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 35, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$immune_group) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 35, units = "mm")


# Feature plot for EDA
add_feat <- "Itgam"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)
