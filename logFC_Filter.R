### calculate ATAC-seq fold change in every bin 
suppressMessages(library(e1071))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(data.table))
suppressMessages(library(xlsx))
suppressMessages(library(gtools))

source("load_BloodData.R")
blood<- loadBlood()
Pan_peaks<- fread("TCGA-ATAC_PanCancer_PeakSet.txt")
colnames(blood)<- Pan_peaks$name
Neutrophil_counts<- blood[1:8 , ]
PBMC_counts<- blood[ 9:259 , ]


set<- set[ , -c(ind_response , ind_id)]
set<- set[ , order(colnames(set))]


source("collapse_samples.R")
Neutrophil_counts<- my_coll(Neutrophil_counts , "Neut")
PBMC_counts<- my_coll(PBMC_counts , "PBMC")


Neut_mean_ratio<- ((colMeans(set[cancerCaseID , ])+1)/(colMeans(Neutrophil_counts)+1) )
PBMC_mean_ratio<- ((colMeans(set[cancerCaseID , ])+1)/(colMeans(PBMC_counts)+1) )

### we assumed that cfDNA drived from Netrophils consist 2/3 of whole plasma cfDNA
blood_mean<-(2*(colMeans(Neutrophil_counts)+1) +(colMeans(PBMC_counts)+1)/3) 
blood_ratio<- ((colMeans(set[cancerCaseID , ])+1) /blood_mean)

colnames(PBMC_counts)<- Pan_peaks$name
colnames(Neutrophil_counts)<- Pan_peaks$name

log2FC_Neut<-log2(Neut_mean_ratio)
log2FC_PBMC<- log2(PBMC_mean_ratio)
log2FC_blood<- log2(blood_ratio)

f1.features<- (names(log2FC_blood))[ which(log2FC_blood >=1.5)]


### get features coordinates
Pan_peaks<- fread("data/TCGA-ATAC_PanCancer_PeakSet.txt")
Pan_peaks$score<-NULL
temp<-Pan_peaks$annotation
Pan_peaks$anno<-temp
Pan_peaks$annotation<-NULL
rm(temp)

Pan_peaks$log2FC_Neutrophil<-log2FC_Neut
Pan_peaks$log2FC_PBMC<-log2FC_PBMC
Pan_peaks$log2FC_blood<- log2FC_blood

source("annotate.R")
f1_coor<- Pan_peaks[Pan_peaks$name%in%f1.features, ]
f1_coor<- my_annotate(f1_coor) ### using annotatePeak adds annotations and gene names 

f1_coor$log2FC_Neutrophil<- log2FC_Neut[f1.features]
f1_coor$log2FC_PBMC<- log2FC_PBMC[f1.features]
f1_coor$mean_ratio_Neut<- Neut_mean_ratio[f1.features]
f1_coor$mean_ratio_PBMC<- PBMC_mean_ratio[f1.features]
f1_coor$mean_ratio_blood<- blood_ratio[f1.features]
f1_coor$log2FC_blood<- log2FC_blood[ f1.features]

write.table(f1_coor , paste(Cancer, "/f1_features.txt" , sep=""), quote = FALSE)
write.csv(f1_coor , file=paste(Cancer , "/f1_features.csv" , sep="")
          , append = FALSE )

#set<- as.data.frame(set)
set <- set[ , f1.features]

