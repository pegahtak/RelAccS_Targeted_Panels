### anotate
suppressMessages(library(ChIPseeker))
suppressMessages(library(plyranges))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
suppressMessages(library(annotate))
library(org.Hs.eg.db)

my_annotate<- function(coor)
{
  
  txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
  G.coor<- as_granges(coor , seqnames=seqnames , start=start , end=end)
  x<- annotatePeak(G.coor , TxDb = txdb ,level = "gene"  , addFlankGeneInfo=TRUE )
  annotations<- as.data.frame(x@anno)
  annotations$symbol<- getSYMBOL(x@anno$geneId , data = 'org.Hs.eg.db')
  annotations$anno<- unlist(map(strsplit(annotations$annotation , split = " "), 1))
  annotations$anno[ grep("exon 1", annotations$annotation)]<-"Exon 1"
  annotations$symbol[which(annotations$anno=="Distal")   ]<- ""
  coor$gene_name<- annotations$symbol
  coor$anno<- annotations$anno
  return(coor)
}
