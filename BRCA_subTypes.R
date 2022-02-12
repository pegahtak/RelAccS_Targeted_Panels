tech_cancerTypes<- unlist( map(strsplit(sequencing_stat$Library_Name , split = "_") ,1))
BRCA_index<- which(tech_cancerTypes=="BRCA")
BRCA_info<- sequencing_stat[ BRCA_index, ]
groups_BRCA<- split(BRCA_info , BRCA_info$BRCA_pam50)

s.names<- unlist(map(strsplit(sample.names , split = "_P") , 1 ))
BRCA_basal_n<- unlist(map(strsplit(unlist(groups_BRCA$Basal[,1]) , split = "_L"), 1))
BRCA_Her2_n<- unlist(map(strsplit(unlist(groups_BRCA$Her2[,1]) , split = "_L"), 1))
BRCA_LumA_n<- unlist(map(strsplit(unlist(groups_BRCA$LumA[,1]) , split = "_L"), 1))
BRCA_LumB_n<- unlist(map(strsplit(unlist(groups_BRCA$LumB[,1]) , split = "_L"), 1))


cancerTypes[which(s.names %in% BRCA_basal_n)]<-"BRCA_Basal"
cancerTypes[which(s.names %in% BRCA_Her2_n)]<-"BRCA_Her2"
cancerTypes[which(s.names %in% BRCA_LumA_n)]<-"BRCA_LumA"
cancerTypes[which(s.names %in% BRCA_LumB_n)]<-"BRCA_LumB"


#otherBRCA_id<-
