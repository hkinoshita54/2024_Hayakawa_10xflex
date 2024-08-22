# Continued from 020_clustering_epi.R
# Without integration
analysis_step <- "050_cell_proportion"

# Load packages ----
library(tidyverse)
library(ggplot2)
library(Seurat)
library(speckle)
library(RColorBrewer)

# Make directories 
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Set palette 
pal <- DiscretePalette(24, palette = "alphabet")

# All ----
fraction <- "all"

## merge 3 fractions
seu_epi <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))
seu_imm <- readRDS(file.path("RDSfiles", "seu_030_imm.RDS"))
seu_str <- readRDS(file.path("RDSfiles", "seu_040_str.RDS"))
seu <- merge(x = seu_epi, y = list(seu_imm, seu_str))
seu$orig.ident <- factor(seu$orig.ident, levels = c("WT", "PT", "PC", "PTC"))
seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Imm.", "Str."))

nclust <- length(levels(seu$cellgroup))
props <- getTransformedProps(clusters = seu$cellgroup, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  # scale_fill_manual(values = pal[1:nclust]) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 2.6, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  # scale_fill_manual(values = pal[1:nclust]) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 2.6, height = 3.2, units = "in", dpi = 150)

# Epithelial ----
fraction <- "epi"
seu <- readRDS(file.path("RDSfiles", "seu_020_epi.RDS"))
nclust <- length(levels(seu$celltype))
props <- getTransformedProps(clusters = seu$celltype, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 2.8, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 2.8, height = 3.2, units = "in", dpi = 150)

# Immune ----
fraction <- "imm"
seu <- readRDS(file.path("RDSfiles", "seu_030_imm.RDS"))
nclust <- length(levels(seu$celltype))
props <- getTransformedProps(clusters = seu$celltype, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 3.7, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 3.7, height = 3.2, units = "in", dpi = 150)

# Stromal ----
fraction <- "str"
seu <- readRDS(file.path("RDSfiles", "seu_040_str.RDS"))
nclust <- length(levels(seu$celltype))
props <- getTransformedProps(clusters = seu$celltype, sample = seu$orig.ident)

prop_df <- as.data.frame(props$Proportions) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(prop_df, file.path(res_path, paste0(fraction, "_prop.xlsx")))

num_df <- as.data.frame(props$Counts) %>% pivot_wider(names_from = sample, values_from = Freq)
openxlsx2::write_xlsx(num_df, file.path(res_path, paste0(fraction, "_num.xlsx")))

ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic()
ggsave(paste0(fraction, "_prop.png"), path = plot_path, width = 3.5, height = 3.2, units = "in", dpi = 150)

ggplot(as.data.frame(props$Counts), aes(x=sample, y=Freq, fill=clusters)) + 
  scale_x_discrete(expand=c(0,0.6), position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal[1:nclust]) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic()
ggsave(paste0(fraction, "_num.png"), path = plot_path, width = 3.5, height = 3.2, units = "in", dpi = 150)
