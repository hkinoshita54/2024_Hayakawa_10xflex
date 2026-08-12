analysis_step <- "051_combine_str"

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
fibro <- readRDS("RDSfiles/seu_042_fibro.RDS")
ec <- readRDS("RDSfiles/seu_041_ec.RDS")
load("RDSfiles/str_celltype1_names.Rdata")

## subset to all the three immune cells
Fib_names <- Cells(fibro)
EC_names <- Cells(ec)
str_names <- c(EC_names, Fib_names, Myo_names, Peri_names, Adipo_names, Glial_names, Neural_names, Meso_names)
seu <- subset(seu_all, cells = str_names)

## Add annotation
### str groups
seu$celltype1 <- "Fibro."
seu$celltype1[EC_names] <- "EC"
seu$celltype1[Myo_names] <- "Myo."
seu$celltype1[Peri_names] <- "Peri."
seu$celltype1[Adipo_names] <- "Adipo."
seu$celltype1[Glial_names] <- "Glial"
seu$celltype1[Neural_names] <- "Neural"
seu$celltype1[Meso_names] <- "Meso."
seu$celltype1 <- factor(seu$celltype1, levels = c("Fibro.", "EC", "Myo.", "Peri.", "Adipo.", "Glial", "Neural", "Meso."))

### annotation within Fibro. and EC
seu$celltype2 <- as.character(seu$celltype1)
seu$celltype2[Fib_names] <- as.character(fibro$celltype2)
seu$celltype2[EC_names] <- as.character(ec$celltype2)
fibro_levels <- levels(fibro$celltype2)
ec_levels <- levels(ec$celltype2)
other_levels <- c("Myo.", "Peri.", "Adipo.", "Glial", "Neural", "Meso.")
all_levels <- c(fibro_levels, ec_levels, other_levels)
seu$celltype2 <- factor(seu$celltype2, levels = all_levels)

# Clustering & plots ----
seu <- cluster(seu, npcs = 50, res = 1)

DimPlot(seu, group.by = "orig.ident", cols = pal$genotype4, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Genotype")
ggsave("sample.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "celltype1", cols = pal$stromal, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "celltype1")
ggsave("celltype1.pdf", path = plot_path, width = 30, height = 30, units = "mm")

Idents(seu) <- "celltype1"
features <- c("Pdgfra","Col1a1",
              "Pecam1","Egfl7",
              "Myh11","Acta2",
              "Rgs5","Notch3",
              "Adipoq","Plin1",
              "Sox10","Plp1",
              "Kit","Etv1",
              "Wt1","Msln")
features <- c("Pdgfra",
              "Pecam1",
              "Myh11",
              "Rgs5",
              "Plin1",
              "Sox10",
              "Etv1",
              "Wt1")
DotPlot(seu, group.by = "celltype1", features = features, dot.scale = 2.3) + 
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
ggsave("dotplot2.pdf", path = plot_path, width = 42, height = 35, units = "mm")

# save seurat object ----
saveRDS(seu, "RDSfiles/seu_051_str_combined.RDS")

# Cell numbers and proportion plot ----
seu <- readRDS("RDSfiles/seu_051_str_combined.RDS")
fraction <- "str"

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
  scale_fill_manual(values = pal$stromal) +
  theme_panel() +
  labs(x = NULL, y = "Cell Proportions") +
  theme(axis.title.y = element_text(angle = 90))
ggsave(paste0(fraction, "_prop.pdf"), path = plot_path, width = 35, height = 35, units = "mm")

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal$stromal) +
  theme_panel() +
  labs(x = NULL, y = "Cell Numbers") +
  theme(axis.title.y = element_text(angle = 90)) & NoLegend()
ggsave(paste0(fraction, "_num.pdf"), path = plot_path, width = 25, height = 35, units = "mm")

