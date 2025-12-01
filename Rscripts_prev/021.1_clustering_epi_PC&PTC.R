# Continued from 020.1_clustering_epi_2.R
# Without integration
analysis_step <- "021.1_clustering_epi_PC&PTC"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(escape)
library(org.Mm.eg.db)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper functions ----
recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
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


# Clustering ----
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
seu <- subset(seu, subset = celltype %in% c("Epi-PC&PTC"))
seu <- recluster(seu, npcs = 50, res = 2)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
markers <- FindAllMarkers(seu, only.pos = TRUE) # c12 & c19 - immune
seu <- subset(seu, idents = c(12,19), invert = TRUE)
seu <- subset(seu, subset = orig.ident %in% c("WT", "PT"), invert = TRUE)
seu$orig.ident <- factor(seu$orig.ident, levels = c("WT", "PT", "PC", "PTC"))

seu <- recluster(seu, npcs = 50, res = 0.6)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", cols = c("#00BFC4", "#C77CFF")) + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "01_epi_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Lrg1", "Cd38", "Msln", "Mmp7", "Cldn4", "Cxcl5", "Pigr")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- "Top2a"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Add celltype annotation and save the Seurat object ----
seu$celltype_2 <- ""
seu$celltype_2[seu$seurat_clusters %in% c(5)] <- "PC1"
seu$celltype_2[seu$seurat_clusters %in% c(4)] <- "PC2"
seu$celltype_2[seu$seurat_clusters %in% c(1,8)] <- "PC3"
seu$celltype_2[seu$seurat_clusters %in% c(6)] <- "PTC1"
seu$celltype_2[seu$seurat_clusters %in% c(3)] <- "PTC2"
seu$celltype_2[seu$seurat_clusters %in% c(0)] <- "PTC3"
seu$celltype_2[seu$seurat_clusters %in% c(2,7)] <- "PTC4"
seu$celltype_2 <- factor(seu$celltype_2, levels = c("PC1", "PC2", "PC3", "PTC1", "PTC2", "PTC3", "PTC4"))
DimPlot(seu, group.by = "celltype_2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
Idents(seu) <- "celltype_2"

# Check markers interactively when necessary
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1000) %>%
  ungroup() -> markers
openxlsx2::write_xlsx(markers, file.path(res_path, "markers.xlsx"))

# dotplot for top 5 markers in each cluster
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DotPlot(seu, group.by = "celltype_2", features = top5$gene) + RotatedAxis()
# wo PTC4
seu_wo4 <- subset(seu, subset = celltype_2 == "PTC4", invert = T)
DotPlot(seu_wo4, group.by = "celltype_2", features = top5$gene[1:30]) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 9, height = 3.6, units = "in", dpi = 150)


VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("vln_QC_.png", path = plot_path, width = 10, height = 4, units = "in", dpi = 150)
## nFeature_RNA low, no actually positive markers > low quality cells??

# cluster profiler ----
# https://zenn.dev/rchiji/books/fdd68b85675c8d/viewer/d4b475
gene_list <- split(markers$gene, markers$cluster)

for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}
gene_list[[7]] <- NULL    # remove PTC4 from the analysis

GOBP <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="BP", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOBP), file.path(res_path, "GOBP.xlsx"))
dotplot(GOBP) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOBP.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)

GOMF <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="MF", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOMF), file.path(res_path, "GOMF.xlsx"))
dotplot(GOMF) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOMF.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)

kegg <- comparegene_listkegg <- compareCluster(geneClusters = gene_list, fun = "enrichKEGG", organism = "mmu")
kegg <- setReadable(x = kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(kegg), file.path(res_path, "kegg.xlsx"))
dotplot(kegg) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("kegg.png", path = plot_path, width = 5, height = 10, units = "in", dpi = 150)

reactome <- compareCluster(geneClusters = gene_list, fun = "enrichPathway", organism = "mouse", readable = TRUE)
openxlsx2::write_xlsx(data.frame(reactome), file.path(res_path, "reactome.xlsx"))
dotplot(reactome) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("reactome.png", path = plot_path, width = 5, height = 10, units = "in", dpi = 150)

wp <- compareCluster(geneClusters = gene_list, fun = "enrichWP", organism = "Mus musculus")
wp <- setReadable(x = wp, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(wp), file.path(res_path, "wp.xlsx"))
dotplot(wp) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("wp.png", path = plot_path, width = 5, height = 10, units = "in", dpi = 150)

# ss enrichment by escape ----
H <- getGeneSets(species = "Mus musculus", library = "H")
names(H) <- gsub("HALLMARK-", "", names(H))
names(H) <- gsub("-", "_", names(H))

# ssGSEA ----
seu <- runEscape(seu, method = "ssGSEA", gene.sets = H,
                 groups = 5000, min.size = 15, new.assay.name = "ssGSEA_H",
                 BPPARAM = SnowParam(workers = 2))
seu <- performNormalization(seu, assay = "ssGSEA_H", gene.sets = H, 
                            scale.factor = seu$nFeature_RNA)

gs_H <- FindAllMarkers(seu, assay = "ssGSEA_H_normalized", min.pct = 0, logfc.threshold = 0)
openxlsx2::write_xlsx(gs_H, file.path(res_path, "ssGSEA_H.xlsx"))

# plot heatmap without PTC4
heatmapEnrichment(seu[,seu$celltype_2 != "PTC4"], group.by = "ident", assay = "ssGSEA_H_normalized", gene.set.use = "all",
                  scale = TRUE, cluster.rows = TRUE, cluster.columns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "right", legend.direction = "vertical")
ggsave("GSEA_H_hm.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)

# save data ----
saveRDS(seu, file = file.path("RDSfiles", "seu_021.1_epi_PC&PTC.RDS"))


# check it on the whole epithelium ----
seu2 <- seu
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
seu$celltype_2 <- as.character(seu$celltype)
seu$celltype_2[colnames(seu2)] <- as.character(seu2$celltype_2)
seu$celltype_2 <- factor(seu$celltype_2, levels = c(levels(seu$celltype)[-10], levels(seu2$celltype_2)))
Idents(seu) <- "celltype_2"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") & NoAxes() & 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype_2_whole_epi.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)


# check regulons from SCENIC analysis in all epithelial cells ----
# auc_mtx is exported from pySCENIC
seu <- seu2
auc_mtx <- read_tsv(file = file.path("int_data", "pyscenic_nig-sc", "epi_2", "auc_mtx.txt")) %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[,colnames(seu)]
seu[["scenicAUC"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "scenicAUC"

# clustering based on auc_mtx (w/o pca)
seu <- RunUMAP(seu, dims = NULL, features = rownames(seu), reduction = NULL, reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, features = rownames(seu), reduction = NULL, dims = NULL, verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)

DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("scenic_cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "celltype_2", label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("scenic_celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident", , cols = c("#00BFC4", "#C77CFF")) + NoAxes()
ggsave("scenic_sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

Idents(seu) <- "celltype_2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("regulon_celltype_2.xlsx")))

# check RNA markers of scenic identified clusters
Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("scenic_cluster_on_Seurat_UMAP.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("markers_scenic_clusters.xlsx")))

# feature plot
fp_path <- file.path("plot", analysis_step, "feature_plot")
fs::dir_create(fp_path)
Idents(seu) <- "celltype"

add_feat <- "Nfkb2(+)"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # min.cutoff = "q25", max.cutoff = "q75"
) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# save data ----
saveRDS(seu, file = file.path("RDSfiles", "seu_021.1_epi_PC&PTC.RDS"))

# edit rss.csv
rss <- read_csv("int_data/pyscenic_nig-sc/epi_2/rss_epi_PC&PTC.csv")
rss[,1] <- sapply(rss[,1], gsub, pattern = "_", replacement = "")
openxlsx2::write_xlsx(rss, file.path(res_path, paste0("rss_epi_PC&PTC.xlsx")))

# edit regulon_z_epi_PC&PTC.csv
regulon_z <- read_csv("int_data/pyscenic_nig-sc/epi_2/regulon_z_epi_PC&PTC.csv")
regulon_z <- regulon_z[-1]
regulon_z$regulon <- paste0(regulon_z$regulon, "(+)")
openxlsx2::write_xlsx(regulon_z, file.path(res_path, paste0("regulon_z_epi_PC&PTC.xlsx")))
