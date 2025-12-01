# Continued from 020.1_clustering_epi_2.R
analysis_step <- "060.1_epi_DE_2"

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
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))

## Dim plot split by sample
Idents(seu) <- "celltype"
DimPlot(seu, label = FALSE, repel = TRUE, cols = "polychrome", split.by = "orig.ident") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

## set genotype_celltype group
seu$gt_ct <- paste(seu$orig.ident, seu$celltype, sep = "_")
gt_ct_lev <- lapply(levels(seu$orig.ident), FUN = paste, levels(seu$celltype), sep = "_") %>% unlist()
seu$gt_ct <- factor(seu$gt_ct, levels = gt_ct_lev)

# DE, Epi-PC&PTC vs Prog., Pre-pit and Neck ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype == "Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PC&PTC vs Progs")
ggsave("PC&PTC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "Epi-PC&PTC"
ref_group <- "Prog_PrePit_Neck"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, PTC vs PC ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct == "PTC_Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$gt_ct == "PC_Epi-PC&PTC"] <- "Ref" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC vs PC")
ggsave("PTC_vs_PC.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)


## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "Epi-PTC"
ref_group <- "Epi-PC"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, Epi-PC vs Prog., Pre-pit and Neck ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct == "PC_Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PC vs Progs")
ggsave("PC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "Epi-PC"
ref_group <- "Prog_PrePit_Neck"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, Epi-PTC vs Prog., Pre-pit and Neck ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct == "PTC_Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC vs Progs")
ggsave("PTC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "Epi-PTC"
ref_group <- "Prog_PrePit_Neck"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, PC&PTC vs WT ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype == "Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$gt_ct %in% c("WT_Prog.", "WT_Pit", "WT_Neck") & seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PC&PTC vs WT")
ggsave("PC&PTC_vs_WT.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PC&PTC"
ref_group <- "WT"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, PC vs WT ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct == "PC_Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$gt_ct %in% c("WT_Prog.", "WT_Pit", "WT_Neck") & seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PC vs WT")
ggsave("PC_vs_WT.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PC"
ref_group <- "WT"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, PTC vs WT ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct == "PTC_Epi-PC&PTC"] <- "Test" 
seu$de_group[seu$gt_ct %in% c("WT_Prog.", "WT_Pit", "WT_Neck") & seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC vs WT")
ggsave("PTC_vs_WT.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PTC"
ref_group <- "WT"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}


# DE, PT vs WT ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct %in% c("PT_Prog.", "PT_Pit", "PT_Neck") & seu$seurat_clusters %in% c(0,2,7,13)] <- "Test"
seu$de_group[seu$gt_ct %in% c("WT_Prog.", "WT_Pit", "WT_Neck") & seu$seurat_clusters %in% c(0,2,7,13)] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PT vs WT")
ggsave("PT_vs_WT.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PT"
ref_group <- "WT"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}
