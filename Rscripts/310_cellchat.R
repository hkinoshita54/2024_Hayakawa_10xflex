# Set up environment ----
# cellchat of PTC sample, using cellchatDB including ECM
analysis_step <- "310_cellchat"

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

saveRDS(cellchat, file.path("RDSfiles", "cellchat_310_PTC.RDS"))

# visualization
interaction_df <- CellChatDB.use$interaction
secreted_signaling <- interaction_df$pathway_name[interaction_df$annotation == "Secreted Signaling"] %>% unique()


groupSize <- as.numeric(table(cellchat@idents))
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# mat <- cellchat@net$weight
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

pathways.show <- c("TGFb", "BMP")

# vertex.receiver = seq(1,5) 
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netVisual_bubble(cellchat, remove.isolate = T, 
                 sources.use = c("TumEC"), targets.use = c("CD4-T", "CD8-T", "Treg",  "gdT", "Prolif.T"),
                 # signaling = pathways.show
                 )
ggsave(paste0("bubble", ".png"), path = plot_path, width = 3.5, height = 6, units = "in", dpi = 150)

# netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(10,11,12,13,15), legend.pos.x = 15)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
