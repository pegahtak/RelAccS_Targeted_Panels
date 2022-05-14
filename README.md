
workflow of designing cancer specific panels based on ATAC-seq data sets. This workflow consists of multiple steps. First relative chromatin accessibility is measured in tumors vs blood cells in logFC.R and bloodFilter.R then differntially accessible regions are identified using DiffBind package (DiffBind.R file). We then implement 3 classification algorthims (SVM, RandomForest and LASSO) on TCGA dataset (classification.R).

To run this code make a folder named data and download TCGA ATAC count matrix from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG then run main.R
