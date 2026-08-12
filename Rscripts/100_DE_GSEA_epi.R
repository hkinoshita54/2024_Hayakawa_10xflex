# Continued from 020_clustering_epi.R
analysis_step <- "100_DE_GSEA_epi"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(msigdbr)
library(fgsea)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# preparation for gsea ----
## gene set
collections <- list()
collections$H <- msigdbr(species = "Mus musculus", category = "H")
collections$KEGG <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
collections$REACTOME <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
collections$BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

## yap signature
yap_list <- readRDS("RDSfiles/yap_list.RDS")
yap_list <- yap_list[4:5]


## Helper function to run fgsea with collapsePathways
run_fgsea <- function(gene_set, gene_set_name){
  fgseaRes <- fgsea(pathways = gene_set, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
  fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], gene_set, rank)
  fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_", gene_set_name, "_", description, ".xlsx")))
}

# Load data ----
seu <- readRDS("RDSfiles/seu_053_epi_combined.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

## Dim plot split by sample
Idents(seu) <- "celltype2"
DimPlot(seu, label = FALSE, repel = TRUE, cols = "polychrome", split.by = "orig.ident") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

## set genotype_celltype group
seu$gt_ct <- paste(seu$orig.ident, seu$celltype2, sep = "_")
gt_ct_lev <- lapply(levels(seu$orig.ident), FUN = paste, levels(seu$celltype2), sep = "_") %>% unlist()
seu$gt_ct <- factor(seu$gt_ct, levels = gt_ct_lev)

# DE, Epi-PTC vs Stem + Prog. ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype2 == "Epi-PTC"] <- "Test" 
seu$de_group[seu$gt_ct %in% c("WT_Stem", "WT_Prog.")] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC vs Progs")
ggsave("PTC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PTC"
ref_group <- "WT_Stem&Prog"
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

### yap signature
fgseaRes <- fgsea(pathways = yap_list, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_yap_", description, ".xlsx")))

### plot result

ggplot(fgseaRes, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = pal$genotype[5]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group)) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_YAP_", description, ".png"), path = plot_path, width = 2.5, height = 0.7 + nrow(fgseaRes)/8, units = "in", dpi = 300) 


# DE, PTC vs PC ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype2 == "Epi-PTC"] <- "Test" 
seu$de_group[seu$celltype2 == "Epi-PC"] <- "Ref" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC vs PC")
ggsave("PTC_vs_PC.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)


## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PTC"
ref_group <- "PC"
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

### yap signature
fgseaRes <- fgsea(pathways = yap_list, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_yap_", description, ".xlsx")))

### plot result

ggplot(fgseaRes, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = pal$genotype[5]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group)) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_YAP_", description, ".png"), path = plot_path, width = 2.5, height = 0.7 + nrow(fgseaRes)/8, units = "in", dpi = 300) 


# DE, PC vs Stem + Prog. ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype2 == "Epi-PC"] <- "Test" 
seu$de_group[seu$gt_ct %in% c("WT_Stem", "WT_Prog.")] <- "Ref" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PC vs Progs")
ggsave("PC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
res <- presto::wilcoxauc(X = seu[["RNA"]]$data, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PC"
ref_group <- "WT_Stem&Prog"
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

### yap signature
fgseaRes <- fgsea(pathways = yap_list, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_yap_", description, ".xlsx")))

### plot result

ggplot(fgseaRes, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = pal$genotype[5]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group)) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_YAP_", description, ".png"), path = plot_path, width = 2.5, height = 0.7 + nrow(fgseaRes)/8, units = "in", dpi = 300) 


# DE, PT vs WT ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$gt_ct %in% c("PT_Stem", "PT_Prog.")] <- "Test"
seu$de_group[seu$gt_ct %in% c("WT_Stem", "WT_Prog.")] <- "Ref"
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


### yap signature
fgseaRes <- fgsea(pathways = yap_list, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_yap_", description, ".xlsx")))

### plot result

ggplot(fgseaRes, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = pal$genotype[5]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group)) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_YAP_", description, ".png"), path = plot_path, width = 2.5, height = 0.7 + nrow(fgseaRes)/8, units = "in", dpi = 300) 
