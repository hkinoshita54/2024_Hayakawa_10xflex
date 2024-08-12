# Last updated 2024-07-16
# Packages for starting scRNA-seq analyses by Seurat

# First set up github following https://happygitwithr.com/github-acct:
# Register at github https://github.com
# Install git with xcode command line tools (run "xcode-select --install" in the shell)

# Set github PAT following https://gist.github.com/Z3tt/3dab3535007acf108391649766409421
# Or, do "git push" in Rstudio, then enter PAT as the password:

# Install the following
options(timeout=3600)
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("devtools")
install.packages('remotes')
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("Matrix")
install.packages("data.table")
install.packages('Seurat')
install.packages("harmony")
install.packages("BiocManager")
install.packages("msigdbr")
BiocManager::install("glmGamPoi")
BiocManager::install("DESeq2")
BiocManager::install("fgsea")
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
devtools::install_github("immunogenomics/presto")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')