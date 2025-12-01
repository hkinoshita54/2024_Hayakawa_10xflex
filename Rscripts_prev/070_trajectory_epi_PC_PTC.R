# Continued from 021_clustering_epi_PC_PTC.R
# Without integration
analysis_step <- "070_trajectory_epi_PC_PTC"

# load packages ----
library(tidyverse)
library(Seurat)
library(CytoTRACE)
library(CytoTRACE2)
library(monocle3)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# load data
seu <- readRDS(file.path("RDSfiles", "seu_021_epiPTC.RDS"))

# CytoTRACE ----
counts_matrix <- LayerData(seu, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu@meta.data$celltype_2)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu@meta.data$celltype_2)
names(pheno) <- colnames(seu)
umap_df <- seu@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = pheno, emb = umap_df, outputDir = res_path)

# CytoTRACE2 ----
seu <- cytotrace2(seu, is_seurat = TRUE, slot_type = "counts", species = 'mouse')
add_feat <- "CytoTRACE2_Relative"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75",
            label = FALSE, repel = TRUE
) + scale_color_viridis_c(option = "B") + NoAxes()
ggsave(paste0(add_feat, ".png"), path = plot_path, width = 5, height = 4, units = "in", dpi = 150)

## the results of CytoTRACE and CytoTRACE2 were in opposite direction...
## maybe better to do it with PTC alone (w/o PC)

################################################################################

# selecting cells by using monocle3 function choose_cells()
cds <- as.cell_data_set(seu)
cds <- choose_cells(cds)
seu <- seu[,cds@colData@rownames]
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "exp_group") + NoAxes()
ggsave("group.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, split.by = "exp_group") + NoAxes()
ggsave("split.png", path = path, width = 12, height = 4, units = "in", dpi = 150)
DimPlot(seu, group.by = "Phase") + NoAxes()
ggsave("cellcycle.png", path = path, width = 5, height = 5, units = "in", dpi = 150)


saveRDS(seu, file = RDSfile)


####################
#### feature plots ####
features <- readLines("gene_set/annotation/01_epi_markers.txt")
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
# check dot plot ----
DotPlot(seu, features = features) + RotatedAxis()
ggsave("dotplot.png", path = path, width = 9, height = 4, units = "in", dpi = 150)
features2 = features[c(2:7,9,10,8,21,23)]
DotPlot(seu, features = features2) + RotatedAxis()
ggsave("dotplot2.png", path = path, width = 6, height = 4, units = "in", dpi = 150)


####################
# monocle3 ----
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
cds <- as.cell_data_set(seu)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu@reductions$umap@cell.embeddings

# learn trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           group_label_size = 5, graph_label_size = 0, cell_size = 1) + NoAxes()
ggsave("monocle3.png", path = paste0(path, "/monocle3"), width = 5, height = 5, units = "in", dpi = 150)

# order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "Isthmus1"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 1) + NoAxes()
ggsave("pseudotime.png", path = paste0(path, "/monocle3"), width = 5, height = 5, units = "in", dpi = 150)


