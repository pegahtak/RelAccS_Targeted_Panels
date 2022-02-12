### load data and collapse samples
library(readxl)

### read TCGA pan caner data set 
### obtained form https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
### download link: https://api.gdc.cancer.gov/data/d6e5358d-491c-4776-8043-b5e49b96706e
data<- readRDS("TCGA_panCancer.rds") 
data<- data[ , -1:-7]

sequencing_stat<- read_excel("TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx" , sheet = 3)
colnames(sequencing_stat)<- sequencing_stat[37, ]
sequencing_stat<- sequencing_stat[ -1:-37, ]

### normalize to total library size
data<- (data/ as.numeric(sequencing_stat$Final_DeDup_reads))
