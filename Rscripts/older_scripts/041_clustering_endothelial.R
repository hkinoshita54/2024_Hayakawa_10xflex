# Continued from 040_clustering_str.R
# Without integration
analysis_step <- "041_clustering_endothelial"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(escape)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper function ----
recluster = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
  return(seu)
}

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 4, height = 4, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Load data ----
seu <- readRDS(file.path("RDSfiles", "seu_040_str.RDS"))
seu <- subset(seu, subset = celltype %in% c("BEC", "LEC", "Peri."))

# Clustering ----
seu <- recluster(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
markers <- FindAllMarkers(seu, only.pos = TRUE) # Check markers > c10 seems immune contamination

seu <- subset(seu, idents = c(10), invert = TRUE)
seu <- recluster(seu, npcs = 50, res = 0.4)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))

features <- readLines(file.path("aux_data", "gene_set", "annotation", "03_str_markers.txt"))
sapply(features, save_fp, seu, fp_path)

features <- readLines(file.path("aux_data", "gene_set", "annotation", "endothelial_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("Eng", "Cd38", "Mctp1", "Kit")
sapply(add_feat, save_fp, seu, fp_path)

# Check markers interactively when necessary
markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Mki67"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 4, height = 4, units = "in", dpi = 150)

# annotation ----
seu$celltype_2 <- ""
seu$celltype_2[seu$seurat_clusters %in% c(5)] <- "ArtEC"
seu$celltype_2[seu$seurat_clusters %in% c(4)] <- "VenEC"
seu$celltype_2[seu$seurat_clusters %in% c(1)] <- "CapEC"
seu$celltype_2[seu$seurat_clusters %in% c(0)] <- "TumEC"
seu$celltype_2[seu$seurat_clusters %in% c(2)] <- "LEC"
seu$celltype_2[seu$seurat_clusters %in% c(3)] <- "Peri."
seu$celltype_2 <- factor(seu$celltype_2, levels = c("ArtEC", "VenEC", "CapEC", "TumEC", "LEC", "Peri."))
DimPlot(seu, group.by = "celltype_2", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype_2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
Idents(seu) <- "celltype_2"

# Check and save markers
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1000) %>%
  ungroup() -> markers
openxlsx2::write_xlsx(markers, file.path(res_path, "markers.xlsx"))

# dotplot
features = c("Sema3g", "Ackr1", "Rgcc", "Mctp1", "Kit", "Prox1", "Rgs5")
DotPlot(seu, group.by = "celltype_2", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5, height = 4, units = "in", dpi = 150)

# save data
saveRDS(seu, file = file.path("RDSfiles", "seu_041_endothelial.RDS"))

# cluster profiler ----
# https://zenn.dev/rchiji/books/fdd68b85675c8d/viewer/d4b475
gene_list <- split(markers$gene, markers$cluster)
for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}
gene_list <- gene_list[1:4]    # analyse BEC

GOBP <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="BP", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOBP), file.path(res_path, "GOBP.xlsx"))
dotplot(GOBP) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOBP.png", path = plot_path, width = 4, height = 6, units = "in", dpi = 150)

GOMF <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="MF", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOMF), file.path(res_path, "GOMF.xlsx"))
dotplot(GOMF) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOMF.png", path = plot_path, width = 4, height = 6, units = "in", dpi = 150)

kegg <- comparegene_listkegg <- compareCluster(geneClusters = gene_list, fun = "enrichKEGG", organism = "mmu")
kegg <- setReadable(x = kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(kegg), file.path(res_path, "kegg.xlsx"))
dotplot(kegg) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("kegg.png", path = plot_path, width = 4, height = 8, units = "in", dpi = 150)

reactome <- compareCluster(geneClusters = gene_list, fun = "enrichPathway", organism = "mouse", readable = TRUE)
openxlsx2::write_xlsx(data.frame(reactome), file.path(res_path, "reactome.xlsx"))
dotplot(reactome) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("reactome.png", path = plot_path, width = 4, height = 6, units = "in", dpi = 150)

wp <- compareCluster(geneClusters = gene_list, fun = "enrichWP", organism = "Mus musculus")
wp <- setReadable(x = wp, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(wp), file.path(res_path, "wp.xlsx"))
dotplot(wp) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("wp.png", path = plot_path, width = 4, height = 6, units = "in", dpi = 150)

 # ss enrichment by escape ----
H <- getGeneSets(species = "Mus musculus", library = "H")
names(H) <- gsub("HALLMARK-", "", names(H))
names(H) <- gsub("-", "_", names(H))

# ssGSEA ----
seu <- subset(seu, subset = celltype_2 %in% c("ArtEC", "VenEC", "CapEC", "TumEC"))
seu <- runEscape(seu, method = "ssGSEA", gene.sets = H,
                 groups = 5000, min.size = 15, new.assay.name = "ssGSEA_H",
                 BPPARAM = SnowParam(workers = 2))
seu <- performNormalization(seu, assay = "ssGSEA_H", gene.sets = H, 
                            scale.factor = seu$nFeature_RNA)

gs_H <- FindAllMarkers(seu, assay = "ssGSEA_H_normalized", min.pct = 0, logfc.threshold = 0)
openxlsx2::write_xlsx(gs_H, file.path(res_path, "ssGSEA_H.xlsx"))

# plot heatmap without PTC4
heatmapEnrichment(seu, group.by = "ident", assay = "ssGSEA_H_normalized", gene.set.use = "all",
                  scale = TRUE, cluster.rows = TRUE, cluster.columns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "right", legend.direction = "vertical")
ggsave("GSEA_H_hm.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)

# check it on the whole epithelium ----
seu2 <- seu
seu <- readRDS(file.path("RDSfiles", "seu_040_str.RDS"))
seu$celltype_2 <- as.character(seu$celltype)
seu$celltype_2[colnames(seu2)] <- as.character(seu2$celltype_2)
seu$celltype_2 <- factor(seu$celltype_2, levels = c(levels(seu$celltype)[-c(9,10,12)], levels(seu2$celltype_2)))
Idents(seu) <- "celltype_2"
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") & NoAxes() & 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype_2_whole_epi.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
