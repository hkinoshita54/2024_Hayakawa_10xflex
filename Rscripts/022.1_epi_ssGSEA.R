# Continued from 020.1_clustering_epi_2.R
# Without integration
analysis_step <- "022.1_epi_ssGSEA"

# load packages ----
library(tidyverse)
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

# Helper function ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

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

# load file and check dimplot ----
seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
Idents(seu) <- "celltype"
DimPlot(seu, group.by = "celltype", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()


# cluster profiler ----
# https://zenn.dev/rchiji/books/fdd68b85675c8d/viewer/d4b475
markers <- FindAllMarkers(seu, only.pos = TRUE)
gene_list <- split(markers$gene, markers$cluster)

for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}
# gene_list[[7]] <- NULL    # remove selected clusters from further analyses

GOBP <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="BP", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOBP), file.path(res_path, "GOBP.xlsx"))
dotplot(GOBP) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOBP.png", path = plot_path, width = 6, height = 8, units = "in", dpi = 150)

GOMF <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont="MF", OrgDb = "org.Mm.eg.db", readable = TRUE)
openxlsx2::write_xlsx(data.frame(GOMF), file.path(res_path, "GOMF.xlsx"))
dotplot(GOMF) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("GOMF.png", path = plot_path, width = 6, height = 8, units = "in", dpi = 150)

kegg <- comparegene_listkegg <- compareCluster(geneClusters = gene_list, fun = "enrichKEGG", organism = "mmu")
kegg <- setReadable(x = kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(kegg), file.path(res_path, "kegg.xlsx"))
dotplot(kegg) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("kegg.png", path = plot_path, width = 6, height = 10, units = "in", dpi = 150)

reactome <- compareCluster(geneClusters = gene_list, fun = "enrichPathway", organism = "mouse", readable = TRUE)
openxlsx2::write_xlsx(data.frame(reactome), file.path(res_path, "reactome.xlsx"))
dotplot(reactome) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("reactome.png", path = plot_path, width = 6, height = 10, units = "in", dpi = 150)

wp <- compareCluster(geneClusters = gene_list, fun = "enrichWP", organism = "Mus musculus")
wp <- setReadable(x = wp, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
openxlsx2::write_xlsx(data.frame(wp), file.path(res_path, "wp.xlsx"))
dotplot(wp) +   theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("wp.png", path = plot_path, width = 6, height = 10, units = "in", dpi = 150)

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

# plot heatmap
heatmapEnrichment(seu[,seu$celltype %in% c("Prog.", "Neck", "Chief", "Pit", "Pariet.", "Epi-PC&PTC")], 
                  group.by = "ident", assay = "ssGSEA_H_normalized", gene.set.use = "all",
                  scale = TRUE, cluster.rows = TRUE, cluster.columns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "right", legend.direction = "vertical")
ggsave("GSEA_H_hm.png", path = plot_path, width = 5, height = 8, units = "in", dpi = 150)

# save data ----
# saveRDS(seu, file = file.path("RDSfiles", "seu_021.1_epi.RDS"))

