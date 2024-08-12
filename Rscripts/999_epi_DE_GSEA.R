# Continued from 020_QC.R
analysis_step <- "030_epi_DE"

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

# Load data ----
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))

# Dim plot
Idents(seu) <- "celltype"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet", split.by = "orig.ident") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 7, height = 3.5, units = "in", dpi = 150)

# DE
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$orig.ident, groups_use = c("T5KA", "T5A_DSS"))
res <- res %>% filter(group == "T5KA" & (pct_in != 0 | pct_out != 0))

test_group <- "T5KA"
ref_group <- "T5A_DSS"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

# GSEA
rank <- res %>% select(feature, auc) %>% deframe()
collections <- list()
collections$H <- msigdbr(species = "Mus musculus", category = "H")
collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

## HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

# C6 with collapsePathways (without plotting)
run_fgsea <- function(gene_set, gene_set_name){
  fgseaRes <- fgsea(pathways = gene_set, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
  fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], gene_set, rank)
  fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_", gene_set_name, "_", description, ".xlsx")))
}

for (i in 2){
  run_fgsea(collections[[i]], names(collections)[i])
}
