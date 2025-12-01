# Continued from 010_QC.R & 020_clustering_epi
# Without integration
analysis_step <- "101_infercnv"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(infercnv)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))
seu <- subset(seu, subset = orig.ident == "WT" & celltype %in% c("Epi-PT", "Epi-PC", "Epi-PTC"), invert = TRUE)
seu <- subset(seu, subset = orig.ident == "PT" & celltype %in% c("Epi-PC", "Epi-PTC"), invert = TRUE)
seu <- subset(seu, subset = orig.ident %in% c("PC", "PTC") & celltype %in% c("Prog.", "Pit", "Neck"), invert = TRUE)
table(seu$celltype, seu$orig.ident)
seu$gt_ct <- paste(seu$orig.ident, seu$celltype, sep = "_")
gt_ct_lev <- lapply(levels(seu$orig.ident), FUN = paste, levels(seu$celltype), sep = "_") %>% unlist()
seu$gt_ct <- factor(seu$gt_ct, levels = gt_ct_lev)

# infercnv ----
cts <- seu[["RNA"]]$counts
anno <- data.frame(cell_name = colnames(seu), celltype = seu$gt_ct)
write_tsv(anno, file = file.path("result", "100_infercnv", "annotation.txt"), col_names = FALSE)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = cts,
                                    annotations_file = file.path("result", "100_infercnv", "annotation.txt"),
                                    delim="\t",
                                    gene_order_file = file.path("aux_data", 
                                                                "mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions"),
                                    ref_group_names = levels(seu$gt_ct)[1:8])
rm(anno, cts)

## start here at night

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = file.path("result", "100_infercnv", "output_dir"),  # dir is auto-created for storing outputs
                             cluster_by_groups = T,   # cluster
                             denoise = T,
                             HMM = T)
seu = infercnv::add_to_seurat(infercnv_output_path = file.path("result", "100_infercnv", "output_dir"),
                                     seurat_obj = seu, # optional
                                     top_n=10)

# Dim plot
DimPlot(seu, group.by = "has_cnv_chr7") + NoAxes() + NoLegend()

# Feature plot
fp_path <- file.path(plot_path, "feature_plot")
add_feat <- ""
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = TRUE, repel = TRUE
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)