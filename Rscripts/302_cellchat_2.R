# cellchat of PTC sample, 
# cellchatDB.new2 (Lrg1-Eng, Grem1/2-Bmpr1a, Rspo1/2/3-Lgr4/5) from 300.2, secreted signaling
# fine clusters
# 2025-12-07
analysis_step <- "302_cellchat_2"

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# helper function ----
cluster = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

# Load data ----
load("RDSfiles/CellChatDB.new2.RData")
# seu <- readRDS("RDSfiles/seu_300_PTC_cellchat.RDS")
seu_all <- readRDS("RDSfiles/seu_052_all_combined.RDS")
Idents(seu_all) <- "celltype2"
seu <- subset(seu_all, subset = orig.ident == "PTC")

# modify group ---
## select cell types to keep
table(seu$celltype2)
keep <- c("Epi-PTC",
          "Treg", "gdT",
          "Conv.-B",
          "Mono", "TAM1", "TAM2", "TAN",
          "iCAF", "myoCAF", "matCAF", "Prolif.CAF",
          "Tum-EC1", "Tum-EC2", "Tum-EC3", "Tum-LEC",
          "Myo.", "Peri.")
seu <- subset(seu, subset = celltype2 %in% keep)
seu$ccgroup <- seu$celltype2
seu$ccgroup <- droplevels(seu$ccgroup)
seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, TAM = c("TAM1", "TAM2"))
# seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, CAF = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"))
seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, `Tum-EC` = c("Tum-EC1", "Tum-EC2", "Tum-EC3"))
table(seu$ccgroup) # at least 50 cells in each category

## check dimplot
seu <- cluster(seu, npcs = 50, res = 1)
Idents(seu) <- "ccgroup"
Idents(seu) %>% levels()
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()  + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("ccgroup.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)
add_feat <- "Acvr1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()

# CellChat analysis ----
## Create cellchat object
# View(CellChatDB.mouse$interaction)
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB.use <- subsetDB(CellChatDB.new2, search = c("Secreted Signaling"), key = "annotation")
cellchat@DB <- CellChatDB.use

## Analyses pipeline following the tutorial
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

df.net <- subsetCommunication(cellchat)
cellchat@LR$LRsig %>% View()

saveRDS(cellchat, file.path("RDSfiles", "cellchat_302_PTC_3.RDS"))

# visualization by dotplot

## Tumor - CAFs
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"))
ggsave(paste0("bubble_PTC_4CAFs", ".png"), path = plot_path, width = 3.5, height = 3.5, units = "in", dpi = 300)

## CAFs - Tumor
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_4CAFs_PTC", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

## test Others are the same as before (i.e. 300_cellchat)
# > same result, save for visualization
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Myo.", "Peri."))
ggsave(paste0("bubble_PTC_Myo", ".png"), path = plot_path, width = 3, height = 3.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Myo.", "Peri."), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_Myo_PTC", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Myo.", "TumEC"), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_TumEC&Myo_PTC", ".png"), path = plot_path, width = 3, height = 5, units = "in", dpi = 300)


######

## CAFs - others
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"), targets.use = c("Treg",  "gdT", "Conv.-B"))
ggsave(paste0("bubble_4CAFs_PTC", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

