# color palette for the project

library(ggplot2)
library(RColorBrewer)
library(khroma)
library(rcartocolor)
library(Polychrome)
library(colorspace)

genotype_col <- c(
  "WT"  = "#999999",
  "PT"  = "#56B4E9",
  "PC"  = "#009E73",
  "TC"  = "#E69F00",
  "PTC" = "#CC79A7"
)

genotype_col4 <- genotype_col[c("WT", "PT", "PC", "PTC")]

cellgroup_col <- c(
  "Epi." = "#D55E00",
  "Imm."     = "#0072B2",
  "Str."    = "#009E73"
)

epi_cols <- brewer.pal(12, "Paired")

tumor_cols <- c(
  # Dys series (cool, PC-related)
  "Dys.1" = "#1B9E77",  # Green-teal (distinct but close to PC anchor)
  "Dys.2" = "#40BFC1",  # Cyan (visually far from Dys.1)
  "Dys.3" = "#5E88C4",  # Blue-leaning teal (cool family, strongly distinct)
  
  # Transitional cluster
  "Trans." = "#7570B3", # Blue-violet, ideal intermediate
  
  # Ca series (warm, PTC-related)
  "Ca.1" = "#F1A2C8",   # Light rose pink
  "Ca.2" = "#E76BF3",   # Violet-magenta (categorically strong)
  "Ca.3" = "#CC79A7",   # PTC genotype color (anchor)
  "Ca.4" = "#A31E70"    # Deep magenta/purple (strong separation)
)

immune_group_col <- c(
  "Tcell" = "#7570B3",
  "Bcell" = "#1B9E77",
  "Mye."  = "#D95F02"
)

t_cols <- brewer.pal(8, "Dark2")

b_cols <- brewer.pal(5, "Set1")

mye_cols <- carto_pal(n = 11, name = "Safe")

stromal_cols <- brewer.pal(8, "Set2")

full_glasbey <- glasbey.colors(30)
ec_cols <- full_glasbey[13:(13 + 8 - 1)] %>% unname()

fibro_sub_cols <- brewer.pal(7, "Set1")

pal <- list(
  genotype      = genotype_col,
  genotype4     = genotype_col[c("WT","PT","PC","PTC")],
  cellgroup     = cellgroup_col,
  epi           = epi_cols,
  tum           = tumor_cols,
  immune_group  = immune_group_col,
  T_cells       = t_cols,
  B_cells       = b_cols,
  myeloid       = mye_cols,
  stromal       = stromal_cols,
  endothelial   = ec_cols,
  fibro_sub     = fibro_sub_cols
)

saveRDS(pal, "RDSfiles/color_palette.RDS")
