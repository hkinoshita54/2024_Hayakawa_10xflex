library(tidyverse)
library(Seurat)
library(biomaRt)

s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# create a table to map mouse gene names to human symbol
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  rename(mouse = external_gene_name, human = hsapiens_homolog_associated_gene_name) %>%
  as_tibble() 
bm[bm == ""] <- NA
bm <- bm %>% na.omit()

# convert human gene names to mouse gene names
s.genes <- bm %>%
  filter(human %in% s.genes) %>%
  pull(mouse) %>%
  unique()

g2m.genes <- bm %>%
  filter(human %in% g2m.genes) %>%
  pull(mouse) %>%
  unique()

# save
save(s.genes, g2m.genes, file = file.path("RDSfiles", "cellcycle_genes_mouse.Rdata"))
