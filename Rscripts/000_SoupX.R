analysis_step <- "000_SoupX"
# following below
# https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

# load packages ----
library(tidyverse)
library(ggplot2)
library(Seurat)
library(SoupX)
library(DropletUtils)


# helper functions ----
cluster = function(seu_obj, npcs = 30, res = 0.8){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- FindNeighbors(seu, dims = 1:npcs)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:npcs)
  return(seu)
}

# Set directories ----
data_path <- "data"
out_path <- "out_data"
sample_names <- c("WT", "PT", "PC", "PTC")
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(data_path, out_path, plot_path, res_path))

# Load data ----
filt.matrix <- Read10X_h5(file.path(data_path, sample_names[1], "sample_filtered_feature_bc_matrix.h5"),use.names = T)
raw.matrix <- Read10X_h5(file.path(data_path, sample_names[1], "sample_raw_feature_bc_matrix.h5"),use.names = T)
raw.matrix <- raw.matrix[rownames(filt.matrix),] # only include features in filt.matrix

# SoupX with Seurat clustering ----
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

## Seurat clustering
seu <- CreateSeuratObject(counts = filt.matrix)
seu <- cluster(seu)
meta <- seu[[]]
umap <- seu@reductions$umap@cell.embeddings

## Add cluster to soup.channel
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)

## estimating contamination
soup.channel  <- autoEstCont(soup.channel)
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)

## write out the adjusted matrix
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
DropletUtils:::write10xCounts(file.path(out_path, paste0(sample_names[1], "_soupX_filt")), adj.matrix)

## store soup.channel to check conatmination fraction later
sc_list <- list()
sc_list[[1]] <- soup.channel

# Sanity checks ----
# following below
# https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
# cntSoggy = rowSums(soup.channel$toc > 0)
# cntStrained = rowSums(adj.matrix > 0)
# mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
# mostZeroed
# tail(sort(rowSums(soup.channel$toc > adj.matrix)/rowSums(soup.channel$toc > 0)), n = 20)
# plotChangeMap(soup.channel, adj.matrix, "Pgc")

# loop over the rest of the sample ----
for (i in 2:length(sample_names)) {
  
  # load data
  filt.matrix <- Read10X_h5(file.path(data_path, sample_names[i], "sample_filtered_feature_bc_matrix.h5"),use.names = T)
  raw.matrix <- Read10X_h5(file.path(data_path, sample_names[i], "sample_raw_feature_bc_matrix.h5"),use.names = T)
  raw.matrix <- raw.matrix[rownames(filt.matrix),] # only include features in filt.matrix
  
  # SoupX with Seurat clustering
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  
  ## Seurat clustering
  seu <- CreateSeuratObject(counts = filt.matrix)
  seu <- cluster(seu)
  meta <- seu[[]]
  umap <- seu@reductions$umap@cell.embeddings
  
  ## Add cluster to soup.channel
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  
  ## estimating contamination
  soup.channel  <- autoEstCont(soup.channel)
  # head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
  
  ## write out the adjusted matrix
  # adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  # DropletUtils:::write10xCounts(file.path(out_path, paste0(sample_names[i], "_soupX_filt")), adj.matrix)
  
  ## store soup.channel to check conatmination fraction later
  sc_list[[i]] <- soup.channel
}

# check contamination fractions
genes_of_interest <- c("Lrg1", "Reg3g", "Reg3b", "Cxcl5", "Saa1", "Saa2")
sc_list[[1]]$soupProfile[genes_of_interest, ] # WT
sc_list[[2]]$soupProfile[genes_of_interest, ] # PT
sc_list[[3]]$soupProfile[genes_of_interest, ] # PC
sc_list[[4]]$soupProfile[genes_of_interest, ] # PTC

head(sc_list[[1]]$soupProfile[order(sc_list[[1]]$soupProfile$est, decreasing = TRUE), ], 30) # WT
head(sc_list[[2]]$soupProfile[order(sc_list[[2]]$soupProfile$est, decreasing = TRUE), ], 30) # PT
head(sc_list[[3]]$soupProfile[order(sc_list[[3]]$soupProfile$est, decreasing = TRUE), ], 30) # PC
head(sc_list[[4]]$soupProfile[order(sc_list[[4]]$soupProfile$est, decreasing = TRUE), ], 30) # PTC

