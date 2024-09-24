# Set up environment ----
analysis_step <- "310_cellcheat"

reticulate::use_condaenv("sc-env")
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))


# Load data ----
seu_all <- readRDS(file.path("RDSfiles", "seu_080_combined.RDS"))
Idents(seu_all) <- "celltype"
seu <- subset(seu_all, subset = orig.ident == "PTC")

# CellChat analysis ----
## Create cellchat object
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB) # use all CellChatDB except for "Non-protein Signaling"
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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_310_Epi-PC&PTC.RDS"))

## need to do visualizations
