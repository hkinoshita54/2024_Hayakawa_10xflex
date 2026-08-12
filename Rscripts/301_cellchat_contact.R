# cellchat of PTC sample, cellchatDB cell-cell contact (Default CellChatDB.mouse)
# 12 clusters: Epi-PTC Treg gdT Conv.-B Mono TAM TAN CAF Tum-EC Tum-LEC Myo. Peri.
# 2025-12-06
analysis_step <- "301_cellchat_contact"

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
seu <- readRDS("RDSfiles/seu_300_PTC_cellchat.RDS")
load("RDSfiles/CellChatDB.new2.RData")
pal <- readRDS("RDSfiles/color_palette.RDS")
pal
cols_ptc <- c("#B15928","#66A61E","#666666","#377EB8","#DDCC77","#332288","#661100","#66C2A5","#F1085C","#720055","#8DA0CB","#E78AC3")
names(cols_ptc) <- levels(seu$ccgroup)

DimPlot(seu, group.by = "ccgroup", cols = cols_ptc, raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Clusters_for_CellChat")
ggsave("ccgroup.pdf", path = plot_path, width = 45, height = 45, units = "mm")

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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_301_PTC_contact.RDS"))


## heatmap ----
### select top 10 pathways manually
### save manually from plot panes
N <- 10
top_pathways <- cellchat@netP$pathways[1:N]

netAnalysis_signalingRole_heatmap(
  cellchat, 
  signaling = top_pathways,
  pattern = "outgoing",
  color.use = cols_ptc,
  width = 4,
  height = 4,
  font.size = 10
)

netAnalysis_signalingRole_heatmap(
  cellchat, 
  signaling = top_pathways,
  pattern = "incoming",
  color.use = cols_ptc,
  width = 4,
  height = 4,
  font.size = 10
)

## circle plot
# netVisual_aggregate(
#   cellchat,
#   signaling = "JAM",
#   layout = "circle",
#   color.use = cols_ptc,
#   # sources.use = sources,
#   # targets.use = targets,
#   label.edge = FALSE
# )

# selected dotplot for figures ----
## This is for legend
# netVisual_bubble(cellchat, remove.isolate = T, 
#                  sources.use = c("Epi-PTC"), targets.use = c("Tum-EC"),
#                  font.size = 7) & 
#   theme(text = element_text(family = "Arial", size = 7),
#         axis.text.x = element_text(angle = 0),
#         legend.title = element_text(size = 6),
#         legend.key.height = rel(0.5))
# ggsave(paste0("bubble_PTC_EC2", ".pdf"), path = plot_path, width = 60, height = 60, units = "mm", device = cairo_pdf) 

## For figures
### Tumor - CAF
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("CAF"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_PTC_CAF", ".pdf"), path = plot_path, width = 35, height = 25, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF"), targets.use = c("Epi-PTC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_CAF_PTC", ".pdf"), path = plot_path, width = 30, height = 30, units = "mm") 

### Tumor - EC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Tum-EC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_PTC_EC", ".pdf"), path = plot_path, width = 35, height = 35, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("Epi-PTC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_EC_PTC", ".pdf"), path = plot_path, width = 30, height = 45, units = "mm") 

### CAF - EC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF"), targets.use = c("Tum-EC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_CAF_EC", ".pdf"), path = plot_path, width = 30, height = 30, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("CAF"),
                 font.size = 7) & NoLegend() &theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_EC_CAF", ".pdf"), path = plot_path, width = 35, height = 25, units = "mm") 
