# cellchat of PTC sample, 
# cellchatDB.new2 (Lrg1-Eng, Grem1/2-Bmpr1a, Rspo1/2/3-Lgr4/5) from 300.2, secreted signaling
# 12 clusters: Epi-PTC Treg gdT Conv.-B Mono TAM TAN *1CAF* Tum-EC Tum-LEC Myo. Peri.
# 2025-12-07
analysis_step <- "303_cellchat_3"

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
  theme_panel() & NoAxes() & labs(title = "Clusters for CellChat") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("ccgroup.pdf", path = plot_path, width = 45, height = 30, units = "mm")

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
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
openxlsx2::write_xlsx(df.net, file = file.path(res_path, "res_cellchat.xlsx"))
cellchat@LR$LRsig %>% View()

saveRDS(cellchat, file.path("RDSfiles", "cellchat_303_PTC_3.RDS"))
cellchat <- readRDS("RDSfiles/cellchat_303_PTC_3.RDS")

# visualization ----

## heatmap
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
netVisual_aggregate(
  cellchat,
  signaling = "TGFb",
  layout = "circle",
  color.use = cols_ptc,
  # sources.use = sources,
  # targets.use = targets,
  label.edge = FALSE
)

# dotplot ----

# PTC to other cells
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Epi-PTC")) & NoLegend()
ggsave(paste0("bubble_PTC_PTC", ".png"), path = plot_path, width = 1.8, height = 2.8, units = "in", dpi = 300) 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Tum-EC", "Tum_LEC")) & NoLegend()
ggsave(paste0("bubble_PTC_EC", ".png"), path = plot_path, width = 1.8, height = 2.8, units = "in", dpi = 300) 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Treg", "gdT", "Conv.-B"),) & NoLegend()
ggsave(paste0("bubble_PTC_T&B", ".png"), path = plot_path, width = 2, height = 2.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Mono", "TAM", "TAN")) & NoLegend()
ggsave(paste0("bubble_PTC_Mye", ".png"), path = plot_path, width = 2, height = 2.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("CAF", "Myo.", "Peri.")) & NoLegend()
ggsave(paste0("bubble_PTC_str", ".png"), path = plot_path, width = 2, height = 3.2, units = "in", dpi = 300)

# other cells to PTC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum-LEC"), targets.use = c("Epi-PTC")) & NoLegend()
ggsave(paste0("bubble_EC_PTC", ".png"), path = plot_path, width = 2, height = 2, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), targets.use = c("Epi-PTC")) & NoLegend()
ggsave(paste0("bubble_Str_PTC", ".png"), path = plot_path, width = 2.2, height = 3, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg",  "gdT", "Conv.-B"), targets.use = c("Epi-PTC")) & NoLegend()
# No significant signaling interactions are inferred!

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Mono", "TAM", "TAN"), targets.use = c("Epi-PTC")) & NoLegend()
ggsave(paste0("bubble_Mye_PTC", ".png"), path = plot_path, width = 2.2, height = 2, units = "in", dpi = 300)

# TumEC,LEC
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum_LEC"), targets.use = c("Tum-EC", "Tum_LEC")) & NoLegend()
ggsave(paste0("bubble_EC_EC", ".png"), path = plot_path, width = 1.2, height = 1.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum_LEC"), targets.use = c("Treg", "gdT", "Conv.-B"))
# > No significant signaling interactions are inferred!

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg", "gdT", "Conv.-B"), targets.use = c("Tum-EC", "Tum_LEC"))
# > No significant signaling interactions are inferred!

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum_LEC"), targets.use = c("Mono", "TAM", "TAN")) & NoLegend()
ggsave(paste0("bubble_EC_Mye", ".png"), path = plot_path, width = 2.2, height = 1.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Mono", "TAM", "TAN"), targets.use = c("Tum-EC", "Tum_LEC")) & NoLegend()
ggsave(paste0("bubble_Mye_EC", ".png"), path = plot_path, width = 1.8, height = 2, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC", "Tum_LEC"), targets.use = c("CAF", "Myo.", "Peri.")) & NoLegend()
ggsave(paste0("bubble_EC_Str", ".png"), path = plot_path, width = 2, height = 2, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), targets.use = c("Tum-EC", "Tum-LEC")) & NoLegend()
ggsave(paste0("bubble_Str_EC", ".png"), path = plot_path, width = 2.5, height = 3.5, units = "in", dpi = 300)

# Str & Immune
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), 
                 targets.use = c("Treg",  "gdT", "Conv.-B")) & NoLegend()
ggsave(paste0("bubble_Str_T&B", ".png"), path = plot_path, width = 2.2, height = 1.6, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg",  "gdT", "Conv.-B"), 
                 targets.use = c("CAF", "Myo.", "Peri.")) & NoLegend()
ggsave(paste0("bubble_T&B_Str", ".png"), path = plot_path, width = 1.5, height = 1.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), 
                 targets.use = c("Mono", "TAM", "TAN")) & NoLegend()
ggsave(paste0("bubble_Str_Mye", ".png"), path = plot_path, width = 3, height = 2.5, units = "in", dpi = 300)

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Mono", "TAM", "TAN"), 
                 targets.use = c("CAF", "Myo.", "Peri.")) & NoLegend()
ggsave(paste0("bubble_Mye_Str", ".png"), path = plot_path, width = 2.5, height = 1.6, units = "in", dpi = 300)

# Among Str
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF", "Myo.", "Peri."), 
                 targets.use = c("CAF", "Myo.", "Peri.")) & NoLegend()
ggsave(paste0("bubble_Str_Str", ".png"), path = plot_path, width = 3.5, height = 3.8, units = "in", dpi = 300)

# Among Immune
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Treg",  "gdT", "Conv.-B", "Mono", "TAM", "TAN"), 
                 targets.use = c("Treg",  "gdT", "Conv.-B", "Mono", "TAM", "TAN")) & NoLegend()
ggsave(paste0("bubble_imm_imm", ".png"), path = plot_path, width = 3, height = 2.2, units = "in", dpi = 300)

# selected dotplot for figures ----
## This is for legend
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Tum-EC"),
                 font.size = 7) & 
  theme(text = element_text(family = "Arial", size = 7),
        axis.text.x = element_text(angle = 0),
        legend.title = element_text(size = 6),
        legend.key.height = rel(0.5))
ggsave(paste0("bubble_PTC_EC2", ".pdf"), path = plot_path, width = 60, height = 60, units = "mm", device = cairo_pdf) 

## For figures
netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Tum-EC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_PTC_EC", ".pdf"), path = plot_path, width = 35, height = 45, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("Epi-PTC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_EC_PTC", ".pdf"), path = plot_path, width = 35, height = 15, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("CAF"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_PTC_CAF", ".pdf"), path = plot_path, width = 30, height = 45, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF"), targets.use = c("Epi-PTC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_CAF_PTC", ".pdf"), path = plot_path, width = 30, height = 45, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("CAF"), targets.use = c("Tum-EC"),
                 font.size = 7) & NoLegend() & theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_CAF_EC", ".pdf"), path = plot_path, width = 35, height = 48, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Tum-EC"), targets.use = c("CAF"),
                 font.size = 7) & NoLegend() &theme(axis.text.x = element_text(angle = 0))
ggsave(paste0("bubble_EC_CAF", ".pdf"), path = plot_path, width = 25, height = 25, units = "mm") 



netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Epi-PTC"), targets.use = c("Mono", "TAM", "TAN"),
                 font.size = 7) & NoLegend() &theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("bubble_PTC_Mye", ".pdf"), path = plot_path, width = 35, height = 45, units = "mm") 

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("Mono", "TAM", "TAN"), targets.use = c("Epi-PTC"),
                 font.size = 7) & NoLegend() &theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("bubble_Mye_PTC", ".pdf"), path = plot_path, width = 45, height = 45, units = "mm") 


# check the distribution of key CCC genes in seurat object ----
features <- c("Lrg1", "Eng", "Grem1", "Wnt5a", "Tgfb1", 
              "Spp1", "Mcam", "Cdh5",
              "Pdgfra", "Itgb4", "Itga6",
              "Plau", "Mdk", "Angptl4", "Lgals9")
DotPlot(seu, group.by = "ccgroup", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -2, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot_ccc_genes.pdf", path = plot_path, width = 75, height = 50, units = "mm")
