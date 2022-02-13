###plots.R
### plots 
#Load necessary packages
suppressMessages( library("gplots"))
suppressMessages( library("devtools"))

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#source("heatmap.3.R")
f4_corr<- read_excel(paste(Cancer , "/f4_features.xlsx" , sep="") , sheet=1)
f4_features<- f4_coor$name

# for heatmap we randomly select 70 samples
PBMC_samples<- sample(nrow(PBMC_counts) , 70) 
f4_ind<- which(Pan_peaks$name%in% f4_features)
heatmap_data_pre<- rbind(n.set[ ,f4_features  ], Neutrophil_counts[ ,f4_ind]
                         , PBMC_counts[PBMC_samples , f4_ind] )
rownames(heatmap_data_pre)<- c(rownames(n.set) , c("N1", "N2", "N3", "N4")
                               , paste("PBMC", PBMC_samples , sep = ""))
heatmap_data<- t(heatmap_data_pre)
data_raw<- heatmap_data
heatmap_data<- log2((heatmap_data+5))
#heatmap_data_final<-normalize.quantiles(heatmap_data)
#rownames(heatmap_data_final)<- rownames(heatmap_data)
#Create  color side bars #columns === cell type
cell_color<- c("burlywood4" ,"coral4" ,"darkorchid1","darkorchid1", "darkorchid1",
               "darkorchid1","darkorchid1","cyan", "darkgrey", "magenta"
               , "darkseagreen2", "darkblue", "darkkhaki", "yellow", "aquamarine", "darkorange"
               , "chartreuse4", "cadetblue4", "brown4", "blueviolet", "deeppink2",
               "gold", "darkslategray3", "pink3", "seagreen2", "khaki", "royalblue1")
cell_color<- as.matrix(cell_color , ncol=1)
cancerGroups<- cbind(cancerGroups ,cell_color)
C.ID<-which(cancerGroups[, 1]==Cancer)
clab<-c()
for( j in 1:nrow(cancerGroups))
  clab<- c(clab , rep(cancerGroups[j, 3] , cancerGroups[j, 2]))

clab<- matrix(c(clab , rep("red" , 74)) , ncol=1)

heat_colnames<-unlist( map(strsplit(colnames(data_raw) , split = "_P"), 1))
basal_id<-which(heat_colnames%in% BRCA_basal_n)
Her2_id<- which(heat_colnames%in%BRCA_Her2_n)
LumA_id<-which(heat_colnames%in% BRCA_LumA_n)
LumB_id<-which(heat_colnames%in% BRCA_LumB_n)
otherBRCA_id<-which(cancerTypes=="BRCA")
#otherBRCA_id<- c(otherBRCA_id , which(heat_colnames%in% unlist(groups_BRCA$Normal[,1 ])))

BRCA_group_clab<- rep("white" , ncol(data_raw))
if(Cancer=="BRCA")
{BRCA_group_clab[basal_id]<- "purple"
BRCA_group_clab[Her2_id]<- "blue"
BRCA_group_clab[LumA_id]<- "green"
BRCA_group_clab[LumB_id]<- "yellow"
BRCA_group_clab[otherBRCA_id]<- "red"
}

main_cancers<- setdiff(cancerGroups[,1] , c("BRCA_Basal", "BRCA_LumA", "BRCA_LumB", "BRCA_Her2"))
col_main_cancers<- cancerGroups[cancerGroups[,1]%in%main_cancers, 3]


if(Cancer=="BRCA")
{
  clab<- cbind(clab , BRCA_group_clab)
  colnames(clab)<- c("Cancer Types" , "BRCA subtypes")
}else
  colnames(clab)<- "Cancer Type"

library(RColorBrewer)
mydist=function(c) {dist(c,method="euclidian")}
distfunc <- function(x) dist(x,method="maximum")
myclust=function(c) {hclust(c,method="average")}
#coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#coul<-brewer.pal(11,"RdBu")
suppressMessages(library(heatmaply))
suppressMessages(library(base64))
coul<-heatmaply::cool_warm(1000)  # , "Cool_warm,(Moreland_2009)"
#x<-log2(data_raw+1)
#qx<-normalize.quantiles(t(x))
#qx<-t(qx)
bk <- c(-100,seq(0,100,by=10))
mycols <- c("red",colorRampPalette(colors = c("white","blue"))(length(bk)-2))

par(mar=c( 1, 1, 1, 1))
pdf(file=paste( Cancer , "/heatmap_main.pdf" , sep="") , width = 30 , height = 20)
heatmap.3(heatmap_data ,col=coul   , key=TRUE , keysize=0.5 , scale = "row"
          , Rowv = FALSE , Colv = TRUE , dendrogram="column"
          ,density.info="none", trace="none",labCol=FALSE , labRow = FALSE  ,cexRow=0.8, ColSideColors =clab
          ,ColSideColorsSize =   0.6  , side.height.fraction = 0.6,  colsep = 1
          , hclustfun=myclust
          , distfun=mydist , margins = c(5, 15)  )

if (Cancer=="BRCA"){
  legend(0.001 , 0.7, legend = c(main_cancers ,"Normal blood") ,title = "Cell Types",
         fill=c(col_main_cancers , "red") , border=FALSE, bty="n", y.intersp = 0.9, cex=1.7
         , horiz=FALSE )
  legend(0.001, 0.85 , title = "BRCA subtypes" , legend = c("basal", "Her2", "LumA" ,"LumB", "unknown"),
         fill=c("purple", "blue", "green", "yellow", "red") , border=FALSE, bty="n", y.intersp =0.9 , cex=1.8
         , x.intersp=0.1  , horiz = FALSE)
}else
  legend(0.001 , 0.8, legend = c(main_cancers ,"Normal blood") ,title = "Cell Types",
         fill=c(col_main_cancers , "red") , border=FALSE, bty="n", y.intersp = 0.9, cex=1.7
         , horiz=FALSE )
# brackets(0.95, 0.9, 0.95, 0.695 , lwd=2 ,ticks=NA, type=4, h=0.5)
dev.off()
graphics.off()



###boxplots
suppressMessages(library(ggplot2))
f4_coor<- f4_coor[ order(f4_coor$RF_imp , decreasing = TRUE), ]
#f4_coor<- f4_coor[ order(f4_coor$log2FC_blood , decreasing = TRUE), ]
#bxplt_n<- min(100 , nrow(f4_coor))
bxplt_n<- min(20 , nrow(f4_coor)) 
names_boxplot_features<- f4_coor$name[1:bxplt_n]
caseIndex<- cancerCaseID
othersIndex<- conrtolID
bloodIndex<- seq(from=411 , to=484 , by=1)
boxplot_counts<- matrix( 0 , nrow = bxplt_n , ncol=ncol(data_raw))
for ( k in 1: nrow(boxplot_counts))
  boxplot_counts[k ,]<- data_raw[  which(rownames(data_raw)==names_boxplot_features[k]) , ]
#  data_raw[1:bxplt_n, ]

cType<-rep("" , ncol(boxplot_counts))
for(j in 1:ncol(boxplot_counts))
{
  if( j %in%caseIndex)
    cType[j]<-Cancer
  else
    if(j %in%othersIndex)
      cType[j]<-"other Cancers"
    else
      if(j %in%bloodIndex)
        cType[j]<-"blood"
}
cellType<-rep(cType , bxplt_n)

peak_name<-c()
for ( j in 1:length(names_boxplot_features))
{
  peak_name<- c(peak_name , rep(names_boxplot_features[j] , ncol(data_raw)))
}
counts<-as.vector(t(boxplot_counts))

boxplot_data<-data.frame(peak_name , cellType , counts)    
boxplot_data$peak_name<- factor(boxplot_data$peak_name , levels = f4_coor$name)
boxplot_data$cellType<- factor(boxplot_data$cellType , levels = c("blood" , Cancer , "other Cancers"))
suppressMessages(library(ggplot2))
par(mar=c(1, 1, 1, 1))
pdf(file=paste(Cancer , "/",Cancer,  ".boxplot.pdf" , sep = "") , width = 100 , height = 70)
boxplt<- ggplot(boxplot_data , aes(x= peak_name , y=counts  , fill=cellType )    , col )+ggtitle(Cancer)+
  theme(plot.title = element_text(face = "bold" , size= 90) , axis.text.x = element_text(angle = 90 , size = 90  , face="bold" ),axis.title.x = element_text(size = 200)
        ,axis.title.y = element_text(size = 130) , legend.key.size = unit(20, "cm") , axis.text.y = element_text( size= 90)
        ,legend.title = element_text( size = 170) ,legend.text = element_text( size = 150 ) )+
  geom_boxplot(lwd = 3 , outlier.size = 5  )+
  coord_cartesian(ylim = c(0, 35) , xlim=c(0,bxplt_n+1))+
  labs(x="peak name" , y=" Normalized ATAC-seq counts")
# geom_text( 10 , 50 , Cancer)
print(boxplt)
dev.off()
graphics.off()



