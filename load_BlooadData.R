loadBlood<- function()
{
  
  Neutrophil_counts1 <- fread("data/Neutrophil/PanCancer.txt" ) ### in Pan_peak order 
  PBMC_counts1<- fread("data/PBMC/PanCancer.txt") ### in Pan_peak order
  
  a<- unlist(map ( strsplit(colnames(PBMC_counts1) , split = "[.]"), 1))
  colnames(PBMC_counts1)<- paste(a ,".bam" , sep="" )
  Pan_peaks<- fread("data/TCGA-ATAC_PanCancer_PeakSet.txt")
  
  # remove excess inf ( chr , start , end , width) as they are provided in Pan_peaks
  Neutrophil_counts1<-Neutrophil_counts1[ , -1:-6]
  PBMC_counts1<-PBMC_counts1[ , -1:-6]
  rownames(Neutrophil_counts1)<- Pan_peaks$name
  rownames(PBMC_counts1)<- Pan_peaks$name
 
  # normalize to library size
  Neut_library_size<- read.table("data/Neutrophil/dedup_reads.txt" , header = T)
  PBMC_library_size<- read.table("data/PBMC/dedup_reads.txt" , header = T)
  Neutrophil_counts1<- (Neutrophil_counts1/as.numeric(Neut_library_size$reads))
  PBMC_counts1<- (PBMC_counts1/as.numeric(PBMC_library_size$reads))
  
  blood_counts1<- cbind(Neutrophil_counts1,PBMC_counts1)
  rownames(blood_counts1)<- rownames(Neutrophil_counts1)

  return(t(blood_counts1))
  
}
