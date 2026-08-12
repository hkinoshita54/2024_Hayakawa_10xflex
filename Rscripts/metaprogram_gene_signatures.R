library(tidyverse)
library(Seurat)

# MP gene sets
ws <- "/Users/hiroto/WORKSPACE"
MP_results <- readRDS(file.path(ws, "2026_human_GC_Atlas/RDSfiles/meta_program_results_2.rds"))
MP_list <- MP_results$MP_list
core9 <- intersect(MP_list$MP_4, MP_list$MP_14)
MP_list[["core9"]] <- core9

# create a table to map mouse gene names to human symbol
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(
  attributes = c(
    "external_gene_name",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_orthology_confidence"
  ),
  mart = mart
) %>%
  distinct() %>%
  rename(
    mouse = external_gene_name,
    human = hsapiens_homolog_associated_gene_name,
    type = hsapiens_homolog_orthology_type,
    confidence = hsapiens_homolog_orthology_confidence
  ) %>%
  filter(
    human != "",
    mouse != "",
    type %in% c("ortholog_one2one", "ortholog_one2many")
  )
saveRDS(bm, "RDSfiles/map_hs_mm_strict.RDS")
bm <- readRDS("RDSfiles/map_hs_mm_strict.RDS")

# convert human gene names to mouse gene names
map_hs_mm <- function(gs){
  gs <- bm %>% 
    filter(human %in% gs) %>% 
    pull(mouse) %>% 
    unique()
}

MP_list_mm <- lapply(MP_list, map_hs_mm)

saveRDS(MP_list_mm, "RDSfiles/MP_list_mm.RDS")
