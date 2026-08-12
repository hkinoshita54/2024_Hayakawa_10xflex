# cellchat of PTC sample, 
# cellchatDB.new2 (Lrg1-Eng, Grem1/2-Bmpr1a, Rspo1/2/3-Lgr4/5) from 300.2, secreted signaling
# most fine clusters ver.2:
# 2026-02-05
analysis_step <- "305_cellchat_5"

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
source("Rscripts/theme_panel.R")

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

## all cells combined 
seu_all <- readRDS("RDSfiles/seu_052_all_combined.RDS")
table(seu_all$celltype2)

## tumor cells
tum <- readRDS("RDSfiles/seu_025_ptc.RDS")
table(tum$celltype3.1)

## add most fine annotations
seu_all$celltype4 <- as.character(seu_all$celltype2)
seu_all$celltype4[Cells(tum)] <- as.character(tum$celltype3.1)
seu_all$celltype4 <- factor(seu_all$celltype4, levels = c(levels(tum$celltype3.1), levels(seu_all$celltype2)))

## > subset PTC
seu <- subset(seu_all, subset = orig.ident == "PTC")

# modify group ---
## select cell types with at least 50 cells
table(seu$celltype4)
keep <- seu[[]] %>% 
  count(celltype4) %>% 
  filter(n > 50) %>% 
  pull(celltype4)
seu <- subset(seu, subset = celltype4 %in% keep)
table(seu$celltype4)

## set celltype3 as the groups for cellchat
seu$ccgroup <- seu$celltype4
seu$ccgroup <- droplevels(seu$ccgroup)
table(seu$ccgroup)

## if you merge groups, run something like below:
# seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, TAM = c("TAM1", "TAM2"))
# seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, CAF = c("iCAF", "myoCAF", "matCAF", "Prolif.CAF"))
# seu$ccgroup <- forcats::fct_collapse(seu$ccgroup, `Tum-EC` = c("Tum-EC1", "Tum-EC2", "Tum-EC3"))

## check dimplot
seu <- cluster(seu, npcs = 50, res = 1)
Idents(seu) <- "ccgroup"
DimPlot(seu, group.by = "ccgroup", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Clusters_for_CellChat")
ggsave("ccgroup.pdf", path = plot_path, width = 60, height = 45, units = "mm")

## save as RDS
saveRDS(seu, "RDSfiles/seu_305_PTC_cellchat_fine.RDS")
seu <- readRDS("RDSfiles/seu_305_PTC_cellchat_fine.RDS")

# CellChat analysis ----
## Create cellchat object
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
load("RDSfiles/CellChatDB.new2.RData")
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
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
openxlsx2::write_xlsx(df.net, file = file.path(res_path, "res_cellchat.xlsx"))
cellchat@LR$LRsig %>% View()

saveRDS(cellchat, file.path("RDSfiles", "cellchat_305_PTC_5.RDS"))
cellchat <- readRDS(file.path("RDSfiles", "cellchat_305_PTC_5.RDS"))

####
# inspect the output excell file, select interaction to show by plot
# this most fine annotation tend to make plots too complicated to interpret

# visualization ----

## grouping
group <- levels(cellchat@idents)
group
gr_tum <- group[1:4]
gr_lym <- group[5:8]
gr_mye <- group[9:13]
gr_caf <- group[14:17]
gr_ec <- group[18:21]
gr_misc <- group[22:24]

## overview
N <- 13
top_pathways <- cellchat@netP$pathways[1:N]

netAnalysis_signalingRole_heatmap(
  cellchat, 
  signaling = top_pathways,
  pattern = "outgoing",
  # color.use = cols_ptc,
  width = 8,
  height = 5,
  font.size = 10
)

netAnalysis_signalingRole_heatmap(
  cellchat, 
  signaling = top_pathways,
  pattern = "incoming",
  # color.use = cols_ptc,
  width = 8,
  height = 5,
  font.size = 10
)

## circle plot
netVisual_aggregate(
  cellchat,
  signaling = "CXCL",
  layout = "circle",
  # sources.use = sources,
  # targets.use = targets,
  label.edge = FALSE
)

## heatmap
netAnalysis_signalingRole_network(
  cellchat, 
  signaling = "CXCL", 
  width = 12, height = 2.5, font.size = 10
)

## bubble plot
netVisual_bubble(
  cellchat, 
  remove.isolate = T, 
  sources.use = gr_tum, 
  targets.use = gr_ec
)
ggsave(paste0("bubble_EC_Mye", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)