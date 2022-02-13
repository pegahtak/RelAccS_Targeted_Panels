suppressMessages(library(e1071))
suppressMessages(library(randomForest))
suppressMessages( library(glmnet))
suppressMessages(library(ROCR))
suppressMessages(library(gbm))
suppressMessages(library(xgboost))

set$response<-Cancer_response
set$id<- IDs
ind_response<- which(colnames(set)=="response")
ind_id<- which(colnames(set)=="id")

set.seed(1234)
folds = createFolds(set$id, k = 5)### to get even distribution use set$id if samples are too few use set$response 
cv = lapply(folds, function(x) {
  training_fold = set[-x, ]
  test_fold =set[x, ]
  training_fold<- training_fold[ , -ind_id]
  test_fold<- test_fold[ , -ind_id]
  ind_response<- which(colnames(training_fold)=="response")
  Y_training_fold<-training_fold[,ind_response]
  Y_test_fold<-test_fold[,ind_response]
  Names_training<- rownames(training_fold)
  Names_test<- rownames(test_fold)
  
  classifierRF<-randomForest(training_fold[,-ind_response], factor(Y_training_fold , c(0,1)),
                             xtest=test_fold[ , -ind_response] , ytest=as.factor(test_fold$response) ,
                             ntree=1000 , importance = TRUE)
  
  prr_RF<- prediction(classifierRF$test$votes[,2], test_fold$response)
  cmRF<- classifierRF$test$confusion
  accuracyRF<- ((cmRF[1,1]+cmRF[2,2])/(cmRF[1,1]+cmRF[1,2]+cmRF[2,1]+cmRF[2,2]) )
  misclassifiedRF <- names(classifierRF$test$predicted[ classifierRF$test$predicted!=Y_test_fold ])
  
  
  classifierSVM<-svm( x=training_fold[ , -ind_response] , y= factor(Y_training_fold , levels = c(0,1)) , kernel="linear"
                      ,probability=TRUE, cost= 0.01 , scale=FALSE)
  res_SVM<- predict(classifierSVM , newdata= test_fold[ , -ind_response] ,probability = TRUE
                    , kernel="linear")
  
  probsSVM<- attr(res_SVM, "probabilities")
  SVM_vars<- coef(classifierSVM)
  SVM_vars<- SVM_vars[-1] ### remove intercept
  misclassifiedSVM <-names(res_SVM[res_SVM[1:length(res_SVM)]!=test_fold$response  ])
  cmSVM<- table(Y_test_fold , res_SVM)
  accuracySVM<- ((cmSVM[1,1]+cmSVM[2,2])/(cmSVM[1,1]+cmSVM[1,2]+cmSVM[2,1]+cmSVM[2+2]))
  
  classifierLASSO<- glmnet(as.matrix(training_fold[, -ind_response]) , as.numeric(Y_training_fold),
                           alpha=1 , family = "binomial"  )  #lambda = cv.lasso$lambda.min
  
 
  LASSOCoefs<- coef(classifierLASSO)[ , 100]
  LASSOCoefs<- LASSOCoefs[-1]# remove intercept
  LASSO.features<- rownames(LASSOCoefs[which(LASSOCoefs!=0)])
  LASSO.prr<- predict(classifierLASSO, newx = as.matrix(test_fold[ , -ind_response]), type="response")
  LASSO.prr<- LASSO.prr[,100 ]
  lasso.pr.class<- LASSO.prr
  cmLASSO<- ModelMetrics::confusionMatrix(Y_test_fold , lasso.pr.class , cutoff = 0.5)
  accuracyLASSO<- ( (cmLASSO[1,1]+cmLASSO[2,2])/(cmLASSO[1,1]+cmLASSO[1,2]+cmLASSO[2,1]+cmLASSO[2,2]))
  misclassifiedLASSO<- names(LASSO.prr[which(lasso.pr.class!=Y_test_fold)])
  
  # Boost.model<- xgboost(data=as.matrix(training_fold[ , -ind_response]) , label = Y_training_fold , n.trees=1000
  #                      , nrounds = 5000 , verbose = 0  )
  # boost.importance<- xgb.importance(model= Boost.model)
  # 
  # GBMModel <- train(response~.,data = training_fold,
  #                   method = "gbm",distribution = "multinomial", verbose = F, n.trees=100)
  
  boost.model=gbm(response~.,data = training_fold,n.trees=500, 
                  distribution = "bernoulli", verbose = F, interaction.depth = 1   ) #interaction.depth=4,cv.folds = 5, shrinkage = 0.01)
  boost.imprtance<-varImp(boost.model , numTrees = 500)
  boost.predict <- predict(boost.model , newdata = test_fold , n.trees = 500 , type="response")
  boost.cm<- ModelMetrics::confusionMatrix( Y_test_fold , boost.predict , cutoff = 0.5  )
  b.p<- boost.predict
  b.p[b.p>=0.5]<-1
  b.p[b.p<0.5]<-0
  MisclassifiedBoost<-Names_test[ which(Y_test_fold!=b.p)]
  accuracy.Boost<- ((boost.cm[1,1]+boost.cm[2,2])/(boost.cm[1,1]+boost.cm[1,2]+boost.cm[2,1]+boost.cm[2,2]))
  
  
  out <- list( accuracyRF, cmRF , classifierRF$importance, prr_RF , classifierRF$test$votes 
               , misclassifiedRF , accuracySVM , cmSVM  , res_SVM , probsSVM , misclassifiedSVM , SVM_vars
               ,accuracyLASSO , cmLASSO , LASSOCoefs , LASSO.prr , misclassifiedLASSO,
               accuracy.Boost ,boost.cm,boost.imprtance , boost.predict,MisclassifiedBoost,
               Y_test_fold)
  names(out)<- c("accuracyRF" , "cmRF" , "importanceRF", "prr_RF" , "votesRF" , "misclassRF"
                 , "accuracySVM" , "cmSVM" , "res_SVM" , "probsSVM", "misclassifiedSVM","SVM_VarCoef",
                 "accuracyLASSO" , "cmLASSO" , "LASSOCoefs" , "LASSO.prr" , "misclassifiedLASSO" ,
                 "accuracyBoost","cmBoost" ,"ImportanceBoost", "PredictBoost","MisclassifiedBoost",
                 "Y_test_fold")
  return(out)
  
})
### report k fold results for each fold and each model  ( AUC , accuracy , misclassified) 
### report over all results for each model ; these results are saved in "Cancer/cv/over_all"

set <- set[ , -c(ind_id , ind_response)]
source("report.R")
acc<-report(Cancer , cv)

### compare models and choose the best one 

accuracyRF <- acc$accuracyRF
accuracySVM<- acc$accuracySVM
accuracyLASSO<- acc$accuracyLASSO
accuracyBoost<- acc$accuracyBoost

#### features suggested by all models 
foldN<- length(cv)
RF.Importance<-(map(cv ,3 ))
RF.imp.mean<- (Reduce('+', RF.Importance)/foldN)
RF.imp.mean<- RF.imp.mean[which(RF.imp.mean[,3]>0)  ,]

Boost.importance<- (map(cv ,20 ))
Boost.imp.mean<- (Reduce('+', Boost.importance)/foldN)
Boost.imp.mean<- data.frame(rownames(Boost.imp.mean) ,Boost.imp.mean )
colnames(Boost.imp.mean)<- c("names" , "Overall")

SVM.coef<- (map(cv ,12 ))
SVM.coef.mean<- (Reduce('+', SVM.coef)/foldN)
SVM.coef.mean<- SVM.coef.mean[ order(abs(SVM.coef.mean) , decreasing = TRUE)] ### absolute value of coefs
LASSO.coef<- (map(cv , 15)) 
LASSO.coef.mean<- rep(0 ,length(LASSO.coef[[1]] ) )
names(LASSO.coef.mean)<- names(LASSO.coef[[1]])

for ( j in 1: length(LASSO.coef[[1]]))
{
  t<-unlist(map(LASSO.coef , j))
  count=0
  for ( k in 1:5)
    if(t[k]!=0)
      count<-count+1
  if (count>=2)
    LASSO.coef.mean[j]<- 1
}

RF.imp.mean<- RF.imp.mean[order(RF.imp.mean[ , 3] , decreasing = TRUE) , ]
RF.features.1 <- (rownames(RF.imp.mean))[ 1:ceiling( 0.1*nrow(RF.imp.mean))]
RF.features<- c()
for ( j in 1:length(RF.features.1))
{
  i.peak<- RF.features.1[j]
  if(mean(n.set[cancerCaseID ,i.peak]) >= ( mean(n.set[conrtolID , i.peak]) ) )
    RF.features<- c(RF.features , i.peak)
}

## Boosting 
Boost.imp.mean<- Boost.imp.mean [order(Boost.imp.mean[,2] , decreasing = TRUE),]
Boost.features.1 <- (rownames(Boost.imp.mean))[ 1:ceiling( 0.1*nrow(Boost.imp.mean))]

Boost.features<- c()
for ( j in 1:length(Boost.features.1))
{
  i.peak<- Boost.features.1[j]
  if(mean(n.set[cancerCaseID ,i.peak]) >= ( mean(n.set[conrtolID , i.peak]) ) )
    Boost.features<- c(Boost.features , i.peak)
}

#SVM features 
SVM.coef.mean<-SVM.coef.mean[ order(abs(SVM.coef.mean)  , decreasing = TRUE )]
SVM.nonZero<- SVM.coef.mean[SVM.coef.mean!=0]
SVM.features.1<- names(SVM.coef.mean) 
SVM.features<-c()
for( j in 1:length(SVM.features.1))
{
  i.peak<- SVM.features.1[j]
  if(mean(n.set[cancerCaseID , i.peak]) >=(mean(n.set[conrtolID,i.peak ])))
    SVM.features<- c(SVM.features , i.peak)
}
SVM.features<- SVM.features[1: ceiling(0.1*length(SVM.features))]

LASSO.features.1<- names(LASSO.coef.mean[LASSO.coef.mean==1])
LASSO.features<-c()
for( j in 1:length(LASSO.features.1))
{
  i.peak<- LASSO.features.1[j]
  if(mean(set[cancerCaseID , i.peak])>=(mean(set[conrtolID , i.peak])) )
    LASSO.features<-c(LASSO.features , i.peak)
}


all_features<-c(RF.features, SVM.features, LASSO.features , Boost.features)
f4_features<- unique(all_features[duplicated(all_features)])
comm_features<-intersect(RF.features, intersect(SVM.features ,LASSO.features)) 

### add RF importance?
t.Importance<-(map(cv ,3 ))
t.imp.mean<- (Reduce('+', t.Importance)/foldN)
RF_imp<- rep(0 , nrow(f3_coor))
for( j in 1:nrow(f3_coor))
{
  RF_imp[j]<- t.imp.mean[which(rownames(t.imp.mean)==f3_coor$name[j]) , 3]
}
f3_coor$RF_imp<-RF_imp
f3_coor$SVM_coef<- SVM.coef.mean
###report features
### sort by q value and RF importance 
f4_coor<- f3_coor[f3_coor$name%in%f4_features , ]
f4_coor$SVM_coefs<-SVM.coef.mean[f4_features]

### add mean value in Cancer 
mean_value<- rep( 0 , nrow(f4_coor))
for( k in 1: nrow(f4_coor))
  mean_value[k]<- mean( n.set[cancerCaseID ,f4_coor$name[k] ])
f4_coor$mean_value<- mean_value
median_value<- rep( 0 , nrow(f4_coor))
for( k in 1: nrow(f4_coor))
  median_value[k]<- median( n.set[cancerCaseID ,f4_coor$name[k] ])
f4_coor$median_value<- median_value
### sort by  logFC ,q value 
f4_coor<- f4_coor[order( f4_coor$RF_imp , decreasing = TRUE),  ]
write.xlsx(f4_coor , file= paste(Cancer, "/f4_features.xlsx", sep="")
           , sheetName = "at least 2 models", append = FALSE )

coor_comm_features<- f3_coor[f3_coor$name%in%comm_features, ]
coor_comm_features$RF_importance<- RF.imp.mean[comm_features, 3]
coor_comm_features$SVM_coefs<-SVM.coef.mean[comm_features]
write.xlsx(coor_comm_features , file=paste(Cancer, "/f4_features.xlsx", sep="")
           , sheetName = "all three models", append = TRUE )

### venn diagram 
RF.id<- which(colnames(n.set)%in%RF.features)
SVM.id<- which(colnames(n.set)%in%SVM.features)
LASSO.id<-  which(colnames(n.set)%in%LASSO.features)
Boost.id<- which(colnames(n.set)%in%Boost.features)

suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
myCol <- brewer.pal(4, "Pastel2")

Venn<- venn.diagram(x=list(RF.id , SVM.id , LASSO.id , Boost.id ) , category.names = c("RF" , "SVM " , "LASSO" , "Boosting" )
                    , filename = paste(Cancer, '/f2_venn_diagramm.png', sep="") ,output=TRUE , height =200 , width=200
                    , lty="blank" , fil=myCol , cex=0.15  ,cat.cex = 0.2 )

### reduce set
set<- set[, f4_features]
### report each model ,in_2 and comm featues 
all_features<- unique(all_features)
coor_all<-f3_coor[f3_coor$name%in%all_features , ]
coor_all$RF_importance<-t.imp.mean[all_features, 3]
coor_all$SVM_coefs<-SVM.coef.mean[all_features]

f4_coor_RF<-coor_all[ coor_all$name%in%RF.features]
f4_coor_SVM<-coor_all[ coor_all$name%in%SVM.features]
f4_coor_LASSO<-coor_all[ coor_all$name%in%LASSO.features]

write.xlsx(f4_coor_RF ,file=paste(Cancer , "/f4_features.xlsx", sep="") 
           , sheetName = "RF" , append = TRUE  )
write.xlsx(f4_coor_SVM ,file=paste(Cancer , "/f4_features.xlsx", sep="")
           , sheetName = "SVM" , append = TRUE )
write.xlsx(f4_coor_LASSO,file=paste(Cancer , "/f4_features.xlsx", sep="") 
           , sheetName = "LASSO" , append = TRUE )


# ### report  decrease in features in each step
Filters<-c("no filter", "logFC" , "blood_Filter" ,"DiffBind" ,"in 2 models")
no_peaks<-c(ncol(n.set) , length(f1.features) , length(f2_features) , length(f3_features)  , length(f4_features))
tab<- data.frame(Filters ,no_peaks)
write.xlsx(tab , file=paste(Cancer , "/peaks_in_each_step.xlsx" , sep=""))
# 


