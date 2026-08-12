# cellchat of PTC sample, cellchatDB secreted signaling (default CellChatDB.mouse + Lrg1-Eng L-R pair)
# 12 clusters: Epi-PTC Treg gdT Conv.-B Mono TAM TAN CAF Tum-EC Tum-LEC Myo. Peri.
# 2025-12-06
analysis_step <- "300_cellchat"

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
seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, CAF = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"))
seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, `Tum-EC` = c("Tum-EC1", "Tum-EC2", "Tum-EC3"))
table(seu$ccgroup) # at least 50 cells in each category

## check dimplot
seu <- cluster(seu, npcs = 50, res = 1)
Idents(seu) <- "ccgroup"
Idents(seu) %>% levels()
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes()  + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("ccgroup.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
add_feat <- "Acvr1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

saveRDS(seu, "RDSfiles/seu_300_PTC_cellchat.RDS")

# Update CellChatDB to include Lrg1-Eng
db <- data.frame(
  ligand       = c("Lrg1",   "Lrg1",    "Lrg1"),
  receptor     = c("Eng",    "Tgfbr1",  "Tgfbr2"),
  pathway_name = c("TGFb"),
  annotation   = c("Secreted Signaling"),
  stringsAsFactors = FALSE
)
db <- db[1,]

CellChatDB.new <- CellChat::updateCellChatDB(
  db             = db,
  merged         = TRUE,
  species_target = "mouse"
)

# CellChat analysis ----
## Create cellchat object
# View(CellChatDB.mouse$interaction)
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB.use <- subsetDB(CellChatDB.new, search = c("Secreted Signaling"), key = "annotation")
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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_300_PTC.RDS"))

# visualization by dotplot

# TumEC to immune, Myo., Peri.
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("Treg", "gdT", "Conv.-B"))
# > No significant signaling interactions are inferred!
# ggsave(paste0("bubble_TumEC_T", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("Mono", "TAM", "TAN"))
ggsave(paste0("bubble_TumEC_Mye", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)


netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("CAF", "Myo.", "Peri."))
ggsave(paste0("bubble_TumEC_Str", ".png"), path = plot_path, width = 2.5, height = 4, units = "in", dpi = 300)

# PTC to other cells
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Tum-EC"))
ggsave(paste0("bubble_PTC_TumEC", ".png"), path = plot_path, width = 3, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Treg", "gdT", "Conv.-B"),)
ggsave(paste0("bubble_PTC_T&B", ".png"), path = plot_path, width = 3, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Mono", "TAM", "TAN"))
ggsave(paste0("bubble_PTC_Mye", ".png"), path = plot_path, width = 3, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("CAF", "Myo.", "Peri."))
ggsave(paste0("bubble_PTC_str", ".png"), path = plot_path, width = 3, height = 4.5, units = "in", dpi = 300)

# other cells to PTC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum-LEC"), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_TumEC_PTC", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_CAFs_PTC", ".png"), path = plot_path, width = 3.5, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Myo.", "TumEC"), targets.use = c("Epi-PTC"))
ggsave(paste0("bubble_TumEC&Myo_PTC", ".png"), path = plot_path, width = 3, height = 5, units = "in", dpi = 300)


# Pericyte and TumEC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Myo.", "Peri."), targets.use = c("Tum-EC"))
ggsave(paste0("bubble_Peri_TumEC", ".png"), path = plot_path, width = 3, height = 4, units = "in", dpi = 300)


# Immune to TumEC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Mono", "TAM", "TAN"), targets.use = c("Tum-EC"))
ggsave(paste0("bubble_Mye_EC", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg",  "gdT", "Conv.-B"), targets.use = c("Tum-EC", "Tum-LEC"))
# No significant signaling interactions are inferred!
# ggsave(paste0("bubble_T&B_EC", ".png"), path = plot_path, width = 4, height = 3, units = "in", dpi = 300)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

