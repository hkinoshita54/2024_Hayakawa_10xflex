library(tidyverse)
library(Seurat)

yap <- openxlsx2::read_xlsx(file = "aux_data/signatures.xlsx", sheet = "YAP", col_names = FALSE, row_names = TRUE)
yap <- yap[,-1]
yap_list <- lapply(seq_len(nrow(yap)), function(i) {
  genes <- unlist(yap[i, ], use.names = FALSE)
  genes <- genes[!is.na(genes) & genes != ""]
  as.character(genes)
})
names(yap_list) <- rownames(yap)
YAP1_UP <- yap_list$YAP1_UP %>% unlist()
YAP1_DN <- yap_list$YAP1_DN %>% unlist()

# create a table to map mouse gene names to human symbol
# library(biomaRt)
# mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
# bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
#   distinct() %>%
#   rename(mouse = external_gene_name, human = hsapiens_homolog_associated_gene_name) %>%
#   as_tibble() 
# bm[bm == ""] <- NA
# bm <- bm %>% na.omit()
# saveRDS(bm, "RDSfiles/map_hs_mm.RDS")
bm <- readRDS("RDSfiles/map_hs_mm.RDS")

# convert human gene names to mouse gene names
YAP1_UP <- bm %>%
  filter(human %in% YAP1_UP) %>%
  pull(mouse) %>%
  unique()

YAP1_DN <- bm %>%
  filter(human %in% YAP1_DN) %>%
  pull(mouse) %>%
  unique()

yap_list$YAP1_UP <- YAP1_UP
yap_list$YAP1_DN <- YAP1_DN

# save
saveRDS(yap_list, file = "RDSfiles/yap_list.RDS")
