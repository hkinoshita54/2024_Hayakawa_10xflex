# Continued from 041_clustering_epi.R
analysis_step <- "101_DE_GSEA_ec"

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
seu <- readRDS("RDSfiles/seu_041_ec.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

## Dim plot split by sample
Idents(seu) <- "celltype2"
# DimPlot(seu, label = FALSE, repel = TRUE, cols = "polychrome", split.by = "orig.ident") + NoAxes()
DimPlot(seu, label = FALSE, repel = TRUE, cols = "polychrome") + NoAxes()
DimPlot(seu, label = FALSE, repel = TRUE, group.by = "orig.ident") + NoAxes()
# ggsave("cluster.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

# ## set genotype_celltype group
# seu$gt_ct <- paste(seu$orig.ident, seu$celltype2, sep = "_")
# gt_ct_lev <- lapply(levels(seu$orig.ident), FUN = paste, levels(seu$celltype2), sep = "_") %>% unlist()
# seu$gt_ct <- factor(seu$gt_ct, levels = gt_ct_lev)

# DE, Tum-EC1,2,3 vs Cap-EC ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype2 %in% c("Tum-EC1","Tum-EC2","Tum-EC3")] <- "Test" 
seu$de_group[seu$celltype2 %in% c("Cap-EC")] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "Tum-EC vs non-Tum-EC")
ggsave("PTC_vs_Progs.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
epi.genes <- c("Atp4a","Atp4b","Pgc","Clps","Chia1","Gkn1","Gkn2","Tff1","Tff2",
               "Muc5ac","Muc6","Mucl3","Psca","Lgals2","Cblif","Pdia2","Lypd8")
mtx <- seu[["RNA"]]$data
mtx <- mtx[!(rownames(mtx) %in% epi.genes),]
res <- presto::wilcoxauc(X = mtx, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "Tum-EC"
ref_group <- "non-Tum-EC"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

tum_ec_up <- res %>% 
  arrange(desc(auc)) %>% 
  slice_head(n = 500) %>% 
  pull(feature)

tum_ec_dn <- res %>% 
  arrange((auc)) %>% 
  slice_head(n = 500) %>% 
  pull(feature)

save(tum_ec_up, tum_ec_dn, file = "tum_ec_up_dn.RData")

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}