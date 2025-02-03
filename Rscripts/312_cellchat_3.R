# Set up environment ----
# cellchat of PTC sample, cellchatDB without ECM, combined PTC1-3 as PTC
analysis_step <- "312_cellchat"

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

# modify idents
keep <- levels(seu$celltype_2)[c(13:17,19:22,24,27:30,33,39:42,46,47,49:51)]
seu <- subset(seu, subset = celltype_2 %in% keep)
seu$celltype_2 <- droplevels(seu$celltype_2)
seu$celltype_2 <- forcats::fct_collapse(seu$celltype_2, PTC = c("PTC1", "PTC2", "PTC3"))
levels(seu$celltype_2)
Idents(seu) <- "celltype_2"
Idents(seu) %>% levels()
DimPlot(seu, cols = "polychrome") + NoAxes()  + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))

# CellChat analysis ----
## Create cellchat object
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"), key = "annotation")
cellchat@DB <- CellChatDB.use

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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_312_PTC_3.RDS"))

# visualization by dotplot

# TumEC to immune cells
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"))
ggsave(paste0("bubble_TumEC_T", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("TAM1", "TAM2", "Mono1",  "Mono2", "TAN"))
ggsave(paste0("bubble_TumEC_Mye", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("myoCAF", "matCAF", "iCAF"))
ggsave(paste0("bubble_TumEC_CAF", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

# PTC to other cells
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC"), targets.use = c("TumEC"))
ggsave(paste0("bubble_PTC_TumEC", ".png"), path = plot_path, width = 3, height = 4.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"),)
ggsave(paste0("bubble_PTC_T", ".png"), path = plot_path, width = 3, height = 3.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC"), targets.use = c("TAM1", "TAM2", "Mono1",  "Mono2", "TAN"),)
ggsave(paste0("bubble_PTC_Mye", ".png"), path = plot_path, width = 3.5, height = 5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("PTC"), targets.use = c("myoCAF", "matCAF", "iCAF"),)
ggsave(paste0("bubble_PTC_CAF", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

# other cells to PTC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("myoCAF", "matCAF", "iCAF"), targets.use = c("PTC"))
ggsave(paste0("bubble_CAFs_PTC", ".png"), path = plot_path, width = 3.5, height = 7.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Myo.", "TumEC"), targets.use = c("PTC"))
ggsave(paste0("bubble_TumEC&Myo_PTC", ".png"), path = plot_path, width = 3, height = 5, units = "in", dpi = 300)


# Pericyte and TumEC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("Peri."))
ggsave(paste0("bubble_TumEC_Peri", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Peri."), targets.use = c("TumEC"))
ggsave(paste0("bubble_Peri_TumEC", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)


# Immune
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TAM1", "TAM2"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"))
ggsave(paste0("bubble_TAM_T", ".png"), path = plot_path, width = 4, height = 4, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg",  "gdT"), targets.use = c("TAM1", "TAM2", "Mono1",  "Mono2", "TAN"))
ggsave(paste0("bubble_T_Mye", ".png"), path = plot_path, width = 4, height = 3, units = "in", dpi = 300)
