# Set up environment ----
# cellchat of PTC sample, using cellchatDB with secreted and cell-cell contact (without ECM)
analysis_step <- "311_cellchat"

reticulate::use_condaenv("sc-env")
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1000 * 1024 ^ 2)

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))


# Load data ----
seu_all <- readRDS(file.path("RDSfiles", "seu_080_combined.RDS"))
Idents(seu_all) <- "celltype_2"
seu <- subset(seu_all, subset = orig.ident == "PTC")

# CellChat analysis ----
## Create cellchat object
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"), key = "annotation")
cellchat@DB <- CellChatDB.use
## Lrg1 is not included in the CellChatDB

## Analyses pipeline following the tutorial
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path("RDSfiles", "cellchat_311_PTC_2.RDS"))

# visualization
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"),)
ggsave(paste0("bubble_TumEC_T", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("TAM1", "TAM2", "Mono1",  "Mono2", "TAN"),)
ggsave(paste0("bubble_TumEC_Mye", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC1", "PTC2", "PTC3"), targets.use = c("TumEC"),)
ggsave(paste0("bubble_PTC123_TumEC", ".png"), path = plot_path, width = 3.5, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("myoCAF", "matCAF", "iCAF"), targets.use = c("PTC1", "PTC2", "PTC3"),)
ggsave(paste0("bubble_CAFs_PTC123", ".png"), path = plot_path, width = 5, height = 8, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC1", "PTC2", "PTC3"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"),)
ggsave(paste0("bubble_TumEC_T", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)