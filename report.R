### report  repeated k fold cv results 
### report cm , accuracy , misclassified , AUC
suppressMessages(library(dplyr))
suppressMessages(library(pROC))
suppressMessages(library(plotROC))
suppressMessages(library(purrr))
suppressMessages(library(ROCR))
report<- function( Cancer , cv ) ### cancer : cases[i] 
{
  n=length(cv)### number of folds
  system(paste0("mkdir " , Cancer,Panel_Type , "/cv" ))
  for ( j in 1: n ){
    system(paste0("mkdir " , Cancer, Panel_Type , "/cv/Fold" , j ))
    
    fold=as.character(j)
    # rep=as.character(rep)
    
    write.table(cv[[j]]$cmRF , paste( Cancer,Panel_Type,"/cv/Fold",j ,"/cmRF.Fold" , fold , ".txt"  , sep="" )
                , quote = FALSE )
    write.table(cv[[j]]$cmSVM , paste( Cancer,Panel_Type,"/cv/Fold",j,"/cmSVM.Fold" , fold , ".txt" , sep="" )
                , quote = FALSE )
    write.table(cv[[j]]$cmLASSO , paste (Cancer,Panel_Type,"/cv/Fold",j,"/cmLASSO.Fold" , fold  ,".txt"  , sep="" )
                , quote = FALSE )
    write.table(cv[[j]]$cmBoost , paste (Cancer,Panel_Type,"/cv/Fold",j,"/cmBoost.Fold" , fold  ,".txt"  , sep="" )
                , quote = FALSE )
    
    write.table(cv[[j]]$misclassRF  , paste( Cancer,Panel_Type,"/cv/Fold", j,"/misclassifiedRF.Fold",fold  , ".txt" , sep="" ))
    write.table(cv[[j]]$misclassifiedSVM  , paste( Cancer,Panel_Type,"/cv/Fold", j,"/misclassifiedSVM.Fold",fold ,  ".txt" , sep="" ))
    write.table(cv[[j]]$misclassifiedLASSO  , paste( Cancer,Panel_Type,"/cv/Fold", j,"/misclassifiedLASSO.Fold",fold  , ".txt" , sep="" ))
    write.table(cv[[j]]$MisclassifiedBoost  , paste( Cancer,Panel_Type,"/cv/Fold", j,"/MisclassifiedBoost.Fold",fold  , ".txt" , sep="" ))
    
    
    suppressMessages(roc.RF<-roc(cv[[j]]$Y_test_fold  , cv[[j]]$votesRF[,2]))
    rf.p2<- performance(cv[[j]]$prr_RF,measure="tpr", x.measure="fpr" )
    fileRF=paste( Cancer,Panel_Type,"/cv/Fold", j,"/AUC.RF.Fold." , fold  , ".pdf" , sep="" )
    pdf(file=fileRF  ,width = 4 , height = 4 )
    plotRF<-plot(rf.p2,main=paste(Cancer , " vs. others",  "  fold" , fold," RF" , sep="")
                 ,col=1,lwd=4, cex.lab=1.5 , cex.main=1.1 , cex.lab=1.1)+
      text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(cv[[j]]$accuracyRF , digits = 4) , sep="")  , cex=1 ,font=2)+
      text(x=0.7 , y=0.8  ,label=paste("AUC= ",round( roc.RF$auc, digits = 4) , sep = "") , cex=1 , font=2 )
    
    print(plotRF)
    dev.off()
    
    ### SVM AUC
    p1.svm<- prediction(cv[[j]]$probsSVM[, 2] , cv[[j]]$Y_test_fold)
    p2.svm<- performance(p1.svm,measure="tpr", x.measure="fpr" )
    suppressMessages(SVM.roc<- roc(cv[[j]]$Y_test_fold ,unlist(p1.svm@predictions) ))
    fileSVM=paste( Cancer,Panel_Type,"/cv/Fold", j,"/AUC.SVM.Fold" , fold , ".pdf" , sep="")
    pdf(file=fileSVM  ,width = 4 , height = 4 )
    plotSVM<- plot(p2.svm,main=paste(Cancer , " vs. others",  "  fold" , fold," SVM " , sep="")
                   ,col=1,lwd=4, cex.lab=1.1 , cex.main=1.1)+
      text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(cv[[j]]$accuracySVM , digits = 4) , sep="")  , cex=1 ,font=2)+
      text(x=0.7 , y=0.8  ,label=paste("AUC= ",round(  pROC::auc(SVM.roc), digits = 4) , sep = "") , cex=1 , font=2 )
    print(plotSVM)
    dev.off()
    
    ### LASSO AUC
    l1<- prediction(cv[[j]]$LASSO.prr , cv[[j]]$Y_test_fold)
    l2<- performance(l1,measure="tpr", x.measure="fpr" )
    suppressMessages(LASSO.roc<- roc(cv[[j]]$Y_test_fold ,unlist(l1@predictions) ))
    fileLASSO=paste( Cancer,Panel_Type,"/cv/Fold", j,"/AUC.LASSO.Fold" , fold  , ".pdf" , sep="")
    pdf(file=fileLASSO  ,width = 4 , height = 4 )
    plotLASSO<-plot(l2,main=paste(Cancer , " vs. others",  "  fold" , fold," LASSO " , sep="")
                    ,col=1,lwd=4, cex.lab=1.1 , cex.main=1.1)+
      text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(cv[[j]]$accuracyLASSO , digits = 4) , sep="")  , cex=1 ,font=2)+
      text(x=0.7 , y=0.8  ,label=paste("AUC= ",round(pROC::auc(LASSO.roc), digits = 4) , sep = "") , cex=1 , font=2 )
    
    print(plotLASSO)
    dev.off()
    
    
    ### Boosting AUC 
    
    
    suppressMessages( Boost.roc<- roc(cv[[j]]$Y_test_fold , cv[[j]]$PredictBoost))
    b.p1<- prediction(cv[[j]]$PredictBoost , cv[[j]]$Y_test_fold )
    b.p2<- performance(b.p1,measure="tpr", x.measure="fpr" )
    fileBoost<-paste( Cancer,Panel_Type,"/cv/Fold", j,"/AUC.Boost.Fold" , fold  , ".pdf" , sep="")
    pdf(file=fileBoost  ,width = 4 , height = 4 )
    plotBoost<-plot(b.p2 ,main=paste(Cancer , " vs. others",  "  fold" , fold," Boosting " , sep="")
                    ,col=1,lwd=4, cex.lab=1.1 , cex.main=1.1)+
      text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(cv[[j]]$accuracyBoost , digits = 4) , sep="")  , cex=1 ,font=2)+
      text(x=0.7 , y=0.8  ,label=paste("AUC= ",round(pROC::auc(Boost.roc), digits = 4) , sep = "") , cex=1 , font=2 )
    
    print(plotBoost)
    dev.off()    
    
  }
  
  
  
  ### over all AUC 
  
  Y<- unlist(map(cv, 23))
  listvotesRF<-(map(cv , 5))
  votesRF<- c( listvotesRF[[1]][,2] , listvotesRF[[2]][,2],listvotesRF[[3]][,2], listvotesRF[[4]][,2]
               ,listvotesRF[[5]][,2])
  listprobsSVM<- (map(cv , 10))
  probsSVM<- c(listprobsSVM[[1]][, 2],listprobsSVM[[2]][, 2] , listprobsSVM[[3]][, 2], listprobsSVM[[4]][, 2]
               , listprobsSVM[[5]][, 2])
  lasso.prr<- unlist(map(cv , 16))
  predictBoost<-unlist(map(cv , 21))
  
  ### accuracy
  cmRF<- Reduce('+' , map(cv ,2))
  cmSVM<- Reduce('+' , map(cv , 8))
  cmLASSO<- Reduce('+' , map(cv , 14))
  cmBoost<- Reduce('+' , map(cv ,19))
  
  accuracyRF<- ((cmRF[1,1]+cmRF[2,2])  /(cmRF[1,1]+cmRF[1,2]+cmRF[2,1]+cmRF[2,2]))
  accuracySVM<- ((cmSVM[1,1]+cmSVM[2,2])  /(cmSVM[1,1]+cmSVM[1,2]+cmSVM[2,1]+cmSVM[2,2]))
  accuracyLASSO<- ((cmLASSO[1,1]+cmLASSO[2,2])  /(cmLASSO[1,1]+cmLASSO[1,2]+cmLASSO[2,1]+cmLASSO[2,2]))
  accuracyBoost<- ((cmBoost[1,1]+cmBoost[2,2])  /(cmBoost[1,1]+cmBoost[1,2]+cmBoost[2,1]+cmBoost[2,2]))
  
  suppressMessages(RF.roc<-roc(Y , votesRF))
  p1.rf<- prediction(votesRF , Y)
  p2.rf<- performance(p1.rf , measure = "tpr" , x.measure = "fpr")
  
  p1.svm<- prediction(probsSVM, Y)
  p2.svm<- performance(p1.svm,measure="tpr", x.measure="fpr" )
  suppressMessages(SVM.roc<- roc(Y ,unlist(p1.svm@predictions) ))
  
  l1<- prediction(lasso.prr , Y)
  l2<- performance(l1,measure="tpr", x.measure="fpr" )
  suppressMessages(LASSO.roc<- roc(Y ,lasso.prr ))
  
  suppressMessages(Boost.roc<- roc(Y , predictBoost))
  
  ### AUC plots 
  system(paste("mkdir " , Cancer, Panel_Type , "/cv/over_all" , sep=""))
  
  #RF
  file1=paste(  Cancer, Panel_Type,"/cv/over_all/","AUC.RF" , ".pdf" , sep="" )
  pdf(file=file1  ,width = 5 , height = 5 )
  plotRF<-plot.roc(RF.roc,main=paste(Cancer , " vs. others",  " RF" , sep="")
                   ,col=1,lwd=4, cex.lab=1.5 , cex.main=1.1 , cex.lab=1.1)
  text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(accuracyRF , digits = 4) , sep="")  , cex=1 ,font=2)
  text(x=0.7 , y=0.8  ,label=paste("AUC= ",round(pROC::auc(RF.roc), digits = 4) , sep = "") , cex=1 , font=2 )
  print(plotRF)
  dev.off()
  
  #Boosting
  file4=paste(  Cancer, Panel_Type,"/cv/over_all/","AUC.Boost" , ".pdf" , sep="" )
  pdf(file=file4  ,width = 5 , height = 5 )
  plotBoost<-plot.roc(Boost.roc,main=paste(Cancer , " vs. others",  " Boosting" , sep="")
                      ,col=1,lwd=4, cex.lab=1.5 , cex.main=1.1 , cex.lab=1.1)
  text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(accuracyBoost , digits = 4) , sep="")  , cex=1 ,font=2)
  text(x=0.7 , y=0.8  ,label=paste("AUC= ",round( pROC::auc(Boost.roc), digits = 4) , sep = "") , cex=1 , font=2 )
  print(plotBoost)
  dev.off()
  
  # SVM
  
  file2=paste( Cancer, Panel_Type,"/cv/over_all/","AUC.SVM" , ".pdf" , sep="" )
  pdf(file=file2  ,width = 5 , height = 5 )
  plotSVM<-plot(SVM.roc,main=paste(Cancer , " vs. others",  " SVM" , sep="")
                ,col=1,lwd=4, cex.lab=1.5 , cex.main=1.1 , cex.lab=1.1 )
  text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(accuracySVM , digits = 4) , sep="")  , cex=1 ,font=2)
  text(x=0.7 , y=0.8  ,label=paste("AUC= ",round( pROC::auc(SVM.roc), digits = 4) , sep = "") , cex=1 , font=2 )
  print(plotSVM)
  dev.off()
  
  ###  LASSO 
  
  file3=paste( Cancer, Panel_Type,"/cv/over_all/","AUC.LASSO" , ".pdf" , sep="" )
  pdf(file=file3  ,width = 5 , height = 5 )
  plotLASSO<-plot(l2,main=paste(Cancer , " vs. others",  " LASSO" , sep="")
                  ,col=1,lwd=4, cex.lab=1.5 , cex.main=1.1 , cex.lab=1.1  
  )+
    text(x=0.7 , y=0.65 , label=paste("accuracy=" , round(accuracyLASSO , digits = 4) , sep="")  , cex=1 ,font=2)+
    text(x=0.7 , y=0.8  ,label=paste("AUC= ",round( pROC::auc(LASSO.roc), digits = 4) , sep = "") , cex=1 , font=2 )
  print(plotLASSO)
  dev.off()
  
  out<-list(accuracyRF , accuracySVM , accuracyLASSO , accuracyBoost)
  names(out)<- c("accuracyRF" , "accuracySVM" , "accuracyLASSO" , "accuracyBoost")
  
  ### sensitivity at 95% specificity
  se.rf<- ci.se(RF.roc , specificities=0.95 , conf.level =0.95 , boot.n=1)
  se.svm<- ci.se(SVM.roc , specificities=0.95 , conf.level =0.95, boot.n=1)
  se.la<- ci.se(roc(Y,lasso.prr) , specificities=0.95 , conf.level =0.95 , boot.n=1)
  sensi<- c(se.rf[2] , se.svm[2], se.la[2])
  sensi<-matrix(sensi , ncol = 3)
  colnames(sensi)<- c("RF" , "SVM" , "LASSO")
  write.table(sensi , paste0(Cancer,Panel_Type, "/sensitivites_overall.txt") , sep = "\t" , quote = F , row.names = F)
  
  return(out)
}
