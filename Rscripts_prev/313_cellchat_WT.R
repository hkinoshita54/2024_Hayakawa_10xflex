# Set up environment ----
# cellchat of WT sample, cellchatDB without ECM
analysis_step <- "313_cellchat_WT"

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
seu <- subset(seu_all, subset = orig.ident == "WT")

# modify idents
keep <- levels(seu$celltype_2)[c(1:5,23,25,26,31,35:37,45,49:51)]
seu <- subset(seu, subset = celltype_2 %in% keep)
seu$celltype_2 <- droplevels(seu$celltype_2)
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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_313_WT.RDS"))

# visualization by dotplot

# EC to immune cells
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CapEC"), targets.use = c("NK", "ILC2", "Mac.",  "cDC"))
ggsave(paste0("bubble_CapEC_Imm", ".png"), path = plot_path, width = 3.5, height = 3.5, units = "in", dpi = 300)

# other cells to Epithelial
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Fib1"), targets.use = keep[1:5])
ggsave(paste0("bubble_Fib1_Epi", ".png"), path = plot_path, width = 3.5, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Fib2"), targets.use = keep[1:5])
ggsave(paste0("bubble_Fib2_Epi", ".png"), path = plot_path, width = 3.5, height = 5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Fib3"), targets.use = keep[1:5])
ggsave(paste0("bubble_Fib3_Epi", ".png"), path = plot_path, width = 3.5, height = 4, units = "in", dpi = 300)

# Pericyte and CapEC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CapEC"), targets.use = c("Peri."))
ggsave(paste0("bubble_CapEC_Peri", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Peri."), targets.use = c("CapEC"))
ggsave(paste0("bubble_Peri_CapEC", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)
