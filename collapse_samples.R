### MyCollapse.R
my_coll <- function (count , cType){
  
  
  for ( k in 1:nrow(count))
  {
    count[k,]<- (count[k ,]/ sum(count[k,]))*1000000
  }
  
  
  
  if (cType=="Neut")
  {
    out<- matrix(0 , nrow=4, ncol= ncol(count))
    for ( ind in 1:4 )
    {
      out[ind,]<- (count[2*ind -1, ]+count[2*ind , ])/2
    }
    return(out)
      
  }
  
 
    if (cType=="PBMC")
    {
      meta<- read.delim("data/PBMC/SraRunTable.txt" , sep = ",")
      meta<- meta[ , c("Run", "submitted_subject_id")]
      ids_pbmc<- unique(meta$submitted_subject_id)
      out<- matrix(0 , ncol = ncol(count) , nrow = length(ids_pbmc))
      for ( ind in 1:length(ids_pbmc))
        {
         s<- meta[ which(meta[,2]==ids_pbmc[ind]) , 1]
         r_ind<-which(rownames(count)%in% s)
         if( length(r_ind)>=1)
        { t<- count[r_ind,]
         if (length(r_ind)>1){
         out[ ind,]<- colMeans(t )}else
           out[ind,]<- t}
        }
       return(out)
    }
  
  
}