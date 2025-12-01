analysis_step <- "010_QC"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(DoubletFinder)

# Make directories ----
fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
## sample_name should be identical to the folder name of each sample
sample_table <- read_excel(file.path("out_data", "sample_table.xlsx")) %>% as.data.frame
sample_name <- sample_table$sample_name
nsample <- length(sample_name)
# exp_group <- sample_table$exp_group
rm(sample_table)

ReadData <- function(sample_name){
  data_dir <- file.path("out_data", paste0(sample_name, "_soupX_filt"))
  cts <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(counts = cts, project = sample_name, min.cells = 3, min.features = 100)
  return(seu)
}
seu_list <-lapply(sample_name, ReadData)

# DoubletFinder ----
source(file.path("Rscripts", "011_DoubletFinder.R"))

## merge the list
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)], add.cell.ids = sample_name)
# rm(seu_list)
seu@meta.data <- seu@meta.data[,c(1:3,7)]
seu$orig.ident <- factor(seu$orig.ident, levels = sample_name)
seu$batch <- "2nd"
seu$batch[seu$orig.ident == "PTC"] <- "1st"
seu$batch <- factor(seu$batch, levels = c("1st", "2nd"))

# QC & filter ----
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_unfil.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

seu <- subset(seu, subset = doublet_finder == "Singlet")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_singlet.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

seu <- subset(seu, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 25)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_filt.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

# Save filtered object
seu <- JoinLayers(seu)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
saveRDS(seu, file.path("RDSfiles", "seu_010_filt.RDS"))