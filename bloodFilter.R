### filter2


neut_pbmc<- read.table("data/blood_cells_Final_PeakSet.txt" , header=TRUE)
HUVEC<- read.table("data/HUVEC_finalPeakSet.txt" , header=TRUE)
HUVEC<-HUVEC[HUVEC$score>=5 , ]
blood_peaks<- rbind(HUVEC[ , 1:3] , neut_pbmc[ , 1:3])
G.f1_coor<- as_granges(f1_coor , seqnames=seqnames , start=start , end=end)
G.blood<- as_granges(blood_peaks , seqnames=chr , start=start , end=end)
f2_coor<- f1_coor[!overlapsAny(G.f1_coor ,G.blood , maxgap=-1), ]
f2_features<- f2_coor$name

set<- set[ , f2_features]
write.csv(f2_coor , file= paste(Cancer, "/f2_features.csv", sep="")
          , append = FALSE )
