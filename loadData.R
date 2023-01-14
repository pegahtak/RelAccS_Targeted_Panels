### load data and collapse samples
library(readxl)
library(purrr)
### read TCGA pan caner data set 
### obtained form https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
### download link: https://api.gdc.cancer.gov/data/d6e5358d-491c-4776-8043-b5e49b96706e
data<- readRDS("data/TCGA-ATAC_PanCan_Raw_Counts.rds")  
data<- data[ , -1:-7]
raw_data<- t(data)

sequencing_stat<- read_excel("TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx" , sheet = 3)
colnames(sequencing_stat)<- sequencing_stat[37, ]
sequencing_stat<- sequencing_stat[ -1:-37, ]

### normalize to total library size
data<- 1000000*(data/ as.numeric(sequencing_stat$Final_DeDup_reads))

### collapse technical replicates
s_names<- map(strsplit(colnames(data) , split = "_L" ),1)
s_names<- unlist(s_names)
samples<- as.data.frame(table(s_names))
n.set <- matrix(0 , nrow = nrow(data) , ncol= nrow(samples))
colnames(n.set)<- samples$s_names
for ( i in 1:nrow(samples))
{
  ind<- which(s_names==samples$s_names[i])
  t<-data[,ind]
  if(length(ind)>1)
    n.set[, i]<-rowMeans(t)
  else
    n.set[ , i]<- t
}
rownames(n.set)<- rownames(data)
rm(data)
n.set<- t(n.set)
