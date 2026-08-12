# Continued from 041_clustering_epi.R
analysis_step <- "102_DE_GSEA_ptc"

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
seu <- readRDS("RDSfiles/seu_025_ptc.RDS")
pal <- readRDS("RDSfiles/color_palette.RDS")

## Check overview
Idents(seu) <- "celltype3.1"
DimPlot(seu, label = FALSE, repel = TRUE, cols = "polychrome") + NoAxes()

# DE, PTC3 vs PTC4 ----
## set groups for DE analysis
seu$de_group <- "Other"
seu$de_group[seu$celltype3.1 %in% c("PTC3")] <- "Test" 
seu$de_group[seu$celltype3.1 %in% c("PTC4")] <- "Ref"
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test", "Other"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC3 vs PTC4")
ggsave("PTC3_vs_PTC4.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
# epi.genes <- c("Atp4a","Atp4b","Pgc","Clps","Chia1","Gkn1","Gkn2","Tff1","Tff2",
#                "Muc5ac","Muc6","Mucl3","Psca","Lgals2","Cblif","Pdia2","Lypd8")
mtx <- seu[["RNA"]]$data
# mtx <- mtx[!(rownames(mtx) %in% epi.genes),]
res <- presto::wilcoxauc(X = mtx, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PTC3"
ref_group <- "PTC4"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

ptc3vs4_up <- res %>% 
  arrange(desc(auc)) %>% 
  slice_head(n = 500) %>% 
  pull(feature)

ptc3vs4_dn <- res %>% 
  arrange((auc)) %>% 
  slice_head(n = 500) %>% 
  pull(feature)

save(ptc3vs4_up, ptc3vs4_dn, file = "tum_ec_up_dn.RData")

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
names(collections$H) <- sub("HALLMARK_", "", names(collections$H))
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}



# DE, PTC3 vs other ----
## set groups for DE analysis
seu$de_group <- "Ref"
seu$de_group[seu$celltype3.1 %in% c("PTC3")] <- "Test" 
seu$de_group <- factor(seu$de_group, levels = c("Ref", "Test"))
DimPlot(seu, group.by = "de_group") & NoAxes() & labs(title = "PTC3 vs other")
ggsave("PTC3_vs_other.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## Wilcoxauc test
# epi.genes <- c("Atp4a","Atp4b","Pgc","Clps","Chia1","Gkn1","Gkn2","Tff1","Tff2",
#                "Muc5ac","Muc6","Mucl3","Psca","Lgals2","Cblif","Pdia2","Lypd8")
mtx <- seu[["RNA"]]$data
# mtx <- mtx[!(rownames(mtx) %in% epi.genes),]
res <- presto::wilcoxauc(X = mtx, y = seu$de_group, groups_use = c("Test", "Ref"))
res <- res %>% filter(group == "Test" & (pct_in != 0 | pct_out != 0))

test_group <- "PTC3"
ref_group <- "other"
description <- paste0(test_group, "_vs_", ref_group)
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

# ptc3vs4_up <- res %>% 
#   arrange(desc(auc)) %>% 
#   slice_head(n = 500) %>% 
#   pull(feature)
# 
# ptc3vs4_dn <- res %>% 
#   arrange((auc)) %>% 
#   slice_head(n = 500) %>% 
#   pull(feature)
# 
# save(ptc3vs4_up, ptc3vs4_dn, file = "tum_ec_up_dn.RData")

## GSEA
rank <- res %>% select(feature, auc) %>% deframe()

### HALLMARK
names(collections$H) <- sub("HALLMARK_", "", names(collections$H))
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

### KEGG, REACTOME, GO-BP and C6 with collapsePathways
for (i in 2:5){run_fgsea(collections[[i]], names(collections)[i])}
