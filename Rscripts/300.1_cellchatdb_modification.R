# modify cellchat db
# 
analysis_step <- "300.1_cellchatdb_modification"

library(tidyverse)
library(CellChat)
library(Seurat)
options(stringsAsFactors = FALSE)

## Make directories
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

# helper ----
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
seu <- readRDS("RDSfiles/seu_300_PTC_cellchat.RDS")
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

CellChatDB <- CellChatDB.mouse

# Inspect structure of the interaction and complex tables ----
# View(CellChatDB$interaction)
# View(CellChatDB$complex)
int <- CellChatDB$interaction
comp <- CellChatDB$complex


# check ligand and receptor expression ----
## Ligands
add_ligands <- c("Rspo1", "Rspo2", "Rspo3",
                 "Grem1", "Grem2")
sapply(add_ligands, save_fp, seu, fp_path) 
## > all the 5 ligands are expressed in CAFs at least slighly

## Lgrs ----
lgrs <- c("Lgr4", "Lgr5", "Lgr6")
sapply(lgrs, save_fp, seu, fp_path) 
## > only Lgr4&5 are in some of the tumor cells

## BMP receptors ----
bmpr_names <- int %>%
  filter(pathway_name=="BMP") %>% 
  pull(receptor) %>% 
  unique()
bmpr_df <- comp[bmpr_names,] 
bmprs <- c(bmpr_df$subunit_1, bmpr_df$subunit_2) %>% unique()
sapply(bmprs, save_fp, seu, fp_path) 
## > all the 6 subunits are moderately expressed in tumor cells

# create data.frame to add to cellchatdb ----
lrg1_eng <- data.frame(
  ligand       = "Lrg1",
  receptor     = "Eng",
  pathway_name = "TGFb",
  annotation   = "Secreted Signaling"
)

rspo_lr <- tidyr::expand_grid(
  ligand = c("Rspo1", "Rspo2", "Rspo3"),
  receptor = c("Lgr4", "Lgr5")
) %>%
  mutate(pathway_name = "WNT", annotation   = "Secreted Signaling")

grem_lr <- tidyr::expand_grid(
  ligand = c("Grem1", "Grem2"),
  receptor = bmpr_names
) %>%
  mutate(pathway_name = "BMP", annotation   = "Secreted Signaling")

newLR <- bind_rows(lrg1_eng, rspo_lr, grem_lr)

# add to cellchatdb ----
CellChatDB.new <- CellChat::updateCellChatDB(
  db             = newLR,
  merged         = TRUE,
  species_target = "mouse"
)

save(CellChatDB.new, file = "RDSfiles/CellChatDB.new.RData")
