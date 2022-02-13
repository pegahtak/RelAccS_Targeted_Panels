### 
suppressMessages(library(DiffBind))
suppressMessages(library(plyranges))
peaks<- as_granges(Pan_peaks[Pan_peaks$name%in%f2_features , 1:3] , seqnames= seqnames , start=start , end=end)

source("load_rawBloodData.R")
blood1<- load_rawBlood()
colnames(blood1)<- Pan_peaks$name
Neutrophil_counts1<- blood1[1:8 ,which(colnames(blood1)%in%f2_features) ]
PBMC_counts1<- blood1[ 9:259 ,which(colnames(blood1)%in%f2_features) ]
PBMC_counts1<- t(PBMC_counts1)
Neutrophil_counts1<- t(Neutrophil_counts1)
rownames(PBMC_counts1)<- f2_features
rownames(Neutrophil_counts1)<- f2_features

### load raw data of cancer cell :
#source("loadRawData.R")

sample.names.2<- rownames(raw_data)
feature.names.2<- colnames(raw_data)

cancerTypes.2=unlist( map(strsplit(sample.names.2 , split = "_") ,1))
cancer_counts <- raw_data[ which(cancerTypes.2==Cancer), ]
if(Cancer=="KIDNEY")
  cancer_counts<- raw_data[ c(which(cancerTypes.2=="KIRC") , which(cancerTypes.2=="KIRP")), ]
if(Cancer=="LUNG")
  cancer_counts<- raw_data[ c(which(cancerTypes.2=="LUAD") , which(cancerTypes.2=="LUSC")), ]

cancer_counts<- cancer_counts[ , order(colnames(cancer_counts))]
cancer_counts<- cancer_counts[ , f2_features]

### build count matrix 
countMatrix<- cbind(t(cancer_counts) , PBMC_counts1 , Neutrophil_counts1)

### build sample IDs 

SampleID<- unlist(map(strsplit(rownames(cancer_counts) , split = "_L") , 1))

PBMCsubject<- list.files(path="data/PBMC/subjects/" , pattern = NULL)
PBMC.data<- data.frame(matrix(ncol=2 ,nrow=0))
colnames(PBMC.data)<- c("ID" , "bam")
for( k in 1:length(PBMCsubject))
{
  s<- as.vector(t(read.table(paste("data/PBMC/subjects/" , PBMCsubject[k] , sep="") , header = FALSE)))
  ID<- c()
  bam<- c()
  for( l in 1:length(s))
  {
    ID<- c(ID , PBMCsubject[k])
    bam<- c(bam , s[l])
  }
  d<- data.frame(ID , bam)
  PBMC.data<- rbind(PBMC.data , d)
}
PBMC.data$bam<- unlist(map(strsplit(PBMC.data$bam , split = "/bam/") , 2))
PBMC.data$bam<- paste( unlist(map( strsplit(PBMC.data$bam , split = ".dedup"), 1))  , ".bam" , sep="")
#ind<- setdiff(PBMC.data$bam , rownames(blood[ 9:259 , ]))
#PBMC.data<- PBMC.data[-(which(PBMC.data$bam==ind)) , ]
### add pBMC & Neut ID
SampleID<- c(SampleID , paste("PBMC_",PBMC.data$ID , sep="") , rep(paste("Neut_" , seq(from=1 , to=4 ) , sep=""), each=2 )  )
Condition<- c( rep("Cancer" , nrow(cancer_counts)) , rep("PBMC" , ncol(PBMC_counts1)) , rep("Neut" , ncol(Neutrophil_counts1)))

###

newDBA <- NULL
system("mkdir DiffBind_counts")
for(  k in 1:ncol(countMatrix))
{
  m<- cbind(Pan_peaks[Pan_peaks$name%in%f2_features,1:3] , countMatrix[ , k])
  #colnames(m)<- c("chromosome" ,"start" , "end" , "counts")
  write.table(m , paste("DiffBind_counts/" ,colnames(countMatrix)[k] , ".txt" , sep="") 
            , row.names = FALSE , sep="\t" , quote=FALSE , col.names = FALSE)
  #newDBA<- dba.peakset(NULL , peaks = peaks , sampID=SampleID[k], 
  # condition=Condition[k], counts=paste("DiffBind_counts/" ,colnames(countMatrix)[k] , ".txt" , sep=""))
}
Counts<- paste("DiffBind_counts/" ,colnames(countMatrix) , ".txt" , sep="")
Replicate<- (as.numeric(duplicated(SampleID))+1)
for( k in 1:length(Replicate)-1){
  if (Replicate[k]==2 && Replicate[k+1]==2){
    Replicate[k+1]<-3
  }
}
# Rep_PBMC<- read.table("data/PBMC/replicates.txt")
# Rep_PBMC<- as.vector(t(Rep_PBMC))
# Rep_PBMC<- as.numeric(Rep_PBMC)
# Rep_Neut<- rep(c(1,2), 4)
# Replicate<- c(Replicate , Rep_PBMC , Rep_Neut)

SampleSheet<- data.frame(SampleID , Condition ,Counts , Replicate)
write.csv(SampleSheet ,"SampleSheet.csv" )
rm(blood1)
rm(Neutrophil_counts1)
rm(PBMC_counts1)
DBA2<- dba(sampleSheet = "SampleSheet.csv")

DBA_res<- dba.normalize(DBA = DBA2  , normalize=DBA_NORM_RLE)
DBA_res<- dba.contrast(DBA_res , reorderMeta=list(Condition="Cancer" ))
DBA_res<- dba.analyze(DBA_res , bBlacklist=FALSE , method = DBA_DESEQ2)
diff<- dba.report(DBA_res ,  contrast = 1,  DataType ="DBA_DATA_GRANGES",bNormalized=TRUE )
diff2<- dba.report(DBA_res ,  contrast = 2,  DataType ="DBA_DATA_GRANGES" ,bNormalized=TRUE )

for( k in 1:dim(diff)[1])
{
  a<- (rownames(diff))[k]
  a<- as.numeric(a)
  diff$Chr[k]<- DBA_res$peaks[[1]]$Chr[a]
}
for( k in 1:dim(diff2)[1])
{
  a<- (rownames(diff2))[k]
  a<- as.numeric(a)
  diff2$Chr[k]<- DBA_res$peaks[[1]]$Chr[a]
}

diff_peaks<- diff[ diff$FDR< 0.05, ]
diff2_peaks<- diff2[ diff2$FDR< 0.05, ]

### PBMC 
G.Pan<- as_granges(Pan_peaks , seqnames = seqnames , start=start , end=end)
diff_peaks$Start<- as.numeric(diff_peaks$Start)
diff_peaks$End<- as.numeric(diff_peaks$End)
G.Diff<- as_granges(diff_peaks , seqnames=Chr , start=Start , end=End)
diff_names<-as.data.frame( join_overlap_left(G.Diff , G.Pan , maxgap = -1))
### Neutrophil
diff2_peaks$Start<- as.numeric(diff2_peaks$Start)
diff2_peaks$End<- as.numeric(diff2_peaks$End)
G.Diff2<- as_granges(diff2_peaks , seqnames=Chr , start=Start , end=End)
diff2_names<-as.data.frame( join_overlap_left(G.Diff2 , G.Pan , maxgap = -1))

system("rm DiffBind_counts/*")

f3_features<-unique( diff_names$name) 
set<- set[ , which(colnames(set) %in%f3_features)]
f3_coor<- f2_coor[ which( f2_coor$name%in%f3_features) , ]

rm(DBA_res)
rm(DBA2)
rm(G.Diff)
rm(G.Diff2)

write.csv(f3_coor , file=paste(Cancer , "/f3_features.csv" , sep="")
          , append = FALSE )
