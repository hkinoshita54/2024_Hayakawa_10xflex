# TF activity inference by decoupleR
# https://bioconductor.org/packages/3.19/bioc/vignettes/decoupleR/inst/doc/tf_sc.html

library(Seurat)
library(decoupleR)
library(tidyverse)
library(pheatmap)

seu <- readRDS(file.path("RDSfiles", "seu_020.1_epi.RDS"))
DimPlot(seu, label = TRUE, cols = "polychrome") & NoAxes()


net <- get_collectri(organism='mouse', split_complexes=FALSE)

# Extract the normalized log-transformed counts
mat <- as.matrix(seu[["RNA"]]$data)

# Run ulm
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Extract ulm and store it in tfsulm in seu
seu[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(seu) <- "tfsulm"

# Scale the data
seu <- ScaleData(seu)
seu[["tfsulm"]]$data <- seu[["tfsulm"]]$scale.data

FeaturePlot(seu, features = "Sox4", cols = c("lightgrey","darkred")) & NoAxes() & NoLegend()


n_tfs <- 25
# Extract activities from object as a long dataframe
df <- t(as.matrix(seu[["tfsulm"]]$data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seu)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 



# Run viper
res_viper <- run_viper(mat=mat, net=net, .source='source', .target='target', .mor='mor', minsize = 5)

# Extract ulm and store it in tfsulm in seu
seu[['tfsvpr']] <- res_viper %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(seu) <- "tfsvpr"

# Scale the data
seu <- ScaleData(seu)
seu[["tfsvpr"]]$data <- seu[["tfsvpr"]]$scale.data

FeaturePlot(seu, features = "Sox4", cols = c("lightgrey","darkred")) & NoAxes() & NoLegend()


n_tfs <- 25
# Extract activities from object as a long dataframe
df <- t(as.matrix(seu[["tfsvpr"]]$data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seu)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
