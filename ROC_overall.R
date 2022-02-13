### ROC overall
library(pROC)
library(roccv)
library(data.table)
ROC_all<- function(cv)
{
  
  res_svm<-list()
  res_RF<- list()
  res_LASSO<- list()
  for( k in 1: length(cv))
  {
    n<- length(cv[[k]]$Y_test_fold)
    
    a<- data.frame(cv[[k]]$probsSVM[ ,2] , as.numeric(cv[[k]]$Y_test_fold) , rep( paste0("Fold" , k) , n) )
    b<- data.frame(cv[[k]]$votesRF[ , 2], cv[[k]]$Y_test_fold , rep( paste0("Fold" , k) , n))
    c<- data.frame(cv[[k]]$LASSO.prr , cv[[k]]$Y_test_fold , rep( paste0("Fold" , k) , n))
    colnames(a)<- c("prediction" ,"response" , "method")
    colnames(b)<- c("prediction" ,"response" , "method")
    colnames(c)<- c("prediction" ,"response" , "method")
    res_svm[[k]]<-a
    res_RF[[k]]<- b
    res_LASSO[[k]]<- c
  }
  res_svm<- rbindlist(res_svm)
  res_RF<- rbindlist(res_RF)
  res_LASSO<- rbindlist(res_LASSO)
  all_svm<- cbind(res_svm[,1:2] , rep("SVM_average" , nrow(res_svm)))
  all_RF<- cbind(res_RF[ , 1:2] , rep("RF_average" , nrow(res_RF)))
  all_LASSO<- cbind(res_LASSO[, 1:2] , rep("LASSO_average" , nrow(res_LASSO)))
  colnames(all_svm)<-  c("prediction" ,"response" , "method")                 
  colnames(all_RF)<-  c("prediction" ,"response" , "method")
  colnames(all_LASSO)<-  c("prediction" ,"response" , "method")
  
  avg_svm<- all_svm
  avg_RF<-all_RF
  avg_LASSO<-all_LASSO
  
  all_svm<- rbind(res_svm , all_svm)
  all_RF<- rbind(res_RF , all_RF)
  all_LASSO<- rbind(res_LASSO , all_LASSO)
  
  #svm
  pdf(paste0(Cancer, Panel_Type, "/svmROC.pdf"))
  rocplot(all_svm,main=paste0("roc plot " ,Cancer , " 5 fold cv" )  )
  dev.off()
  graphics.off()
  
  #RF
  pdf(paste0(Cancer, Panel_Type, "/rfROC.pdf"))
  rocplot(all_RF,main=paste0("roc plot " ,Cancer , " 5 fold cv" )  )
  dev.off()
  graphics.off()
  
  #LASSO
  pdf(paste0(Cancer, Panel_Type, "/lassoROC.pdf"))
  rocplot(all_LASSO,main=paste0("roc plot " ,Cancer , " 5 fold cv" )  )
  dev.off()
  graphics.off()
  
  #total
  total<- rbind(avg_svm , avg_RF , avg_LASSO)
  pdf(paste0(Cancer, Panel_Type, "/total_ROC.pdf"))
  rocplot(total,main=paste0("roc plot " ,Cancer , " 5 fold cv" )  )
  dev.off()
  graphics.off()
  
  
  return()
}