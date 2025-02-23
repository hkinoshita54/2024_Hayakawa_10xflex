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
# remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE) # did not work
devtools::install_github("immunogenomics/presto")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# Additional packages
BiocManager::install("escape")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install("sva")
devtools::install_local("/Users/kinoshitahiroto/Downloads/CytoTRACE_0.3.3.tar.gz")
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
devtools::install_github("saeyslab/nichenetr")

devtools::install_github("jinworks/CellChat")
devtools::install_github("renozao/NMF@devel")
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

BiocManager::install("biomaRt")
BiocManager::install("densvis")
devtools::install_github("califano-lab/PISCES")

BiocManager::install("decoupleR")
BiocManager::install("OmnipathR")
BiocManager::install("viper")
