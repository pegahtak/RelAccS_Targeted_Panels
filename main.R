### f1:logfC f2:DEseq f3:HUVEC f4:in 2 models
## main 
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))
suppressMessages(library(data.table))
suppressMessages(library(readr))
source("loadData.R")

sample.names<- rownames(n.set)
feature.names<- colnames(n.set)

cancerTypes=unlist( map(strsplit(sample.names , split = "_") ,1))
set=n.set
set<- as.data.frame(n.set)
source("BRCA_subTypes.R") ## gets subTypes of breast cancer 

### size of cancer types 
cancerGroups<- as.data.frame(table(cancerTypes))
#write.table(cancerGroups ,"cancerTypesInfo.txt" , sep="\t" , row.names = FALSE
#            , quote = FALSE)

IDs<- rep(0 , nrow(set))
for ( j in 1:length(cancerGroups))
{
  IDs[which(cancerTypes==cancerGroups[j])]<-j
}

cases<- c("BRCA", "COAD", "KIDNEY", "LIHC" ,"LUNG" , "PRAD", "STAD" )

Panel_Type="Open"
###open panels
for (  i in 1:length(cases))
{
 
  set=as.data.frame(n.set)
  ### choose cancer Type ( choose i )
  Cancer=cases[i]
  Cancer_response<-rep(0 , nrow(set))
  
  print(paste("*********", Cancer, "*******"))
  if(cases[i]=="KIDNEY")
  {Cancer_response[which(cancerTypes== "KIRP")] <- 1
  Cancer_response[which(cancerTypes== "KIRC")] <- 1
  } else
    
    if(cases[i]=="LUNG")
    {
      Cancer_response[which(cancerTypes== "LUAD")] <- 1
      Cancer_response[which(cancerTypes== "LUSC")] <- 1
    }else
      if( cases[i]=="BRCA")
      {
        
        #Cancer_response[which(cancerTypes=="BRCA_Basal")]<-1
        Cancer_response[which(cancerTypes=="BRCA_LumB")]<-1
        Cancer_response[which(cancerTypes=="BRCA_LumA")]<-1
        Cancer_response[which(cancerTypes=="BRCA_Her2")]<-1
        Cancer_response[which(cancerTypes=="BRCA")]<-1
      }else
        Cancer_response[which(cancerTypes== cases[i])] <- 1
  
  
  
  cancerCaseID<- which(Cancer_response==1)
  conrtolID<- setdiff(seq(from=1 , to=nrow(set) , by=1 ) ,cancerCaseID )
  
  
  set$response<- Cancer_response
  set$response<- factor(set$response , levels = c(0 , 1))
  ind_response<- which(colnames(set)=="response")
  set$id<- IDs
  ind_id<- which(colnames(set)=="id")
  
  system(paste( "mkdir " , Cancer , sep=""))
  system(paste0( "mkdir " , Cancer , Panel_Type))
  source("logFC_Filter.R")
  source("bloodFilter.R")
  source("DiffBind.R")  
  source("classification.R")
  source("ROC_overall.R")
  ROC_all(cv)
  source("plots.R")
 
}

### closed Panels
Panel_Type="Closed"
for (  i in 1:length(cases))
{
  
  set=as.data.frame(n.set)
  ### choose cancer Type ( choose i )
  Cancer=cases[i]
  Cancer_response<-rep(0 , nrow(set))
  
  print(paste("*********", Cancer, "*******"))
  if(cases[i]=="KIDNEY")
  {Cancer_response[which(cancerTypes== "KIRP")] <- 1
  Cancer_response[which(cancerTypes== "KIRC")] <- 1
  } else
    
    if(cases[i]=="LUNG")
    {
      Cancer_response[which(cancerTypes== "LUAD")] <- 1
      Cancer_response[which(cancerTypes== "LUSC")] <- 1
    }else
      if( cases[i]=="BRCA")
      {
        
        Cancer_response[which(cancerTypes=="BRCA_Basal")]<-1
        Cancer_response[which(cancerTypes=="BRCA_LumB")]<-1
        Cancer_response[which(cancerTypes=="BRCA_LumA")]<-1
        Cancer_response[which(cancerTypes=="BRCA_Her2")]<-1
        Cancer_response[which(cancerTypes=="BRCA")]<-1
      }else
        Cancer_response[which(cancerTypes== cases[i])] <- 1
  
  
  
  cancerCaseID<- which(Cancer_response==1)
  conrtolID<- setdiff(seq(from=1 , to=nrow(set) , by=1 ) ,cancerCaseID )
  
  
  set$response<- Cancer_response
  set$response<- factor(set$response , levels = c(0 , 1))
  ind_response<- which(colnames(set)=="response")
  set$id<- IDs
  ind_id<- which(colnames(set)=="id")
  
  
  system(paste0( "mkdir " , Cancer , Panel_Type))
  source("logFC_Filter.R")
  source("bloodFilter.R")
  source("DiffBind.R")  
  source("classification.R")
  source("ROC_overall.R")
  ROC_all(cv)
  source("plots.R")
  
}