# Continued from 041_clustering_endothelial.R
analysis_step <- "062.2_endothelial_DE_GSEA_3"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(msigdbr)
library(fgsea)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# preparation for gsea ----
## gene set
collections <- list()
collections$H <- msigdbr(species = "Mus musculus", category = "H")
collections$KEGG <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
collections$REACTOME <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
collections$BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

## Helper function to run fgsea with collapsePathways
run_fgsea <- function(gene_set, gene_set_name){
  fgseaRes <- fgsea(pathways = gene_set, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
  fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], gene_set, rank)
  fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_", gene_set_name, "_", description, ".xlsx")))
}

# Load data ----
seu <- readRDS(file.path("RDSfiles", "seu_041.1_endothelial.RDS"))

## Dim plot split by sample
Idents(seu) <- "celltype_2"
DimPlot(seu, label = FALSE, repel = TRUE, cols = "alphabet", split.by = "orig.ident") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

## set genotype_celltype group
seu$gt_ct <- paste(seu$orig.ident, seu$celltype_2, sep = "_")
gt_ct_lev <- lapply(levels(seu$orig.ident), FUN = paste, levels(seu$celltype_2), sep = "_") %>% unlist()
seu$gt_ct <- factor(seu$gt_ct, levels = gt_ct_lev)

# DE----
## set groups for DE analysis
seu$de_group <- ""
seu$de_group[seu$celltype_2 %in% c("TumEC")] <- "test" 
seu$de_group[seu$celltype_2 %in% c("CapEC")] <- "ref"
seu$de_group <- factor(seu$de_group, levels = c("ref", "test"))

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("test", "ref"))
res <- res %>% filter(group == "test" & (pct_in != 0 | pct_out != 0))

test_group <- "TumEC"
ref_group <- "CapEC"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% dplyr::select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}



