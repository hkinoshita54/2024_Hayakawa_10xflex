# cellchat of PTC sample, 
# contact signaling
# most fine clusters
# 2026-03-04
analysis_step <- "306_cellchat_6_finer_contact"

library(CellChat)
library(patchwork)
library(Seurat)
source("Rscripts/theme_panel.R")
options(stringsAsFactors = FALSE)

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# Load data ----
seu <- readRDS("RDSfiles/seu_305_PTC_cellchat_fine.RDS")

# CellChat analysis ----
## Create cellchat object
cellchat <- createCellChat(object = seu, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"), key = "annotation")
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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_306_PTC_6.RDS"))

# visualization ----
# inspect the output excell file, select interaction to show by plot
# this most fine annotation tend to make plots too complicated to interpret
cellchat <- readRDS(file.path("RDSfiles", "cellchat_306_PTC_6.RDS"))

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
# netVisual_aggregate(
#   cellchat,
#   signaling = "CXCL",
#   layout = "circle",
#   # sources.use = sources,
#   # targets.use = targets,
#   label.edge = FALSE
# )
# 
# ## heatmap
# netAnalysis_signalingRole_network(
#   cellchat, 
#   signaling = "CXCL", 
#   width = 12, height = 2.5, font.size = 10
# )
# 
# ## bubble plot
# netVisual_bubble(
#   cellchat, 
#   remove.isolate = T, 
#   sources.use = gr_tum, 
#   targets.use = gr_ec
# )
# ggsave(paste0("bubble_EC_Mye", ".png"), path = plot_path, width = 3, height = 3, units = "in", dpi = 300)