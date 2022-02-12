### load data and collapse samples
library(readxl)

### read TCGA pan caner data set 
### obtained form https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
### download link: https://api.gdc.cancer.gov/data/d6e5358d-491c-4776-8043-b5e49b96706e
data<- readRDS("TCGA_panCancer.rds") 
