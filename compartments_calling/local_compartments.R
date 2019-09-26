
library(HiTC)


##Inputs
dirWT = "matrix/" #Hi-C pro output matrix in the matrix format for HiTC 

H3K27ac <- "H3K27ac_keep1.idr_0.01.filtered"

binsize=10000 #binsize of HiC data
type="Iovino_WT"
nbin=as.integer(400/binsize*6400)  ##number of bins to use to calculate local compartments (2.56 Mb)
##
H3K27ac_data = read.table(H3K27ac)
H3K27ac_gr = makeGRangesFromDataFrame(H3K27ac_data,seqnames.field="V1",start.field="V2",end.field="V3",strand.field="V4")



all_files = list.files(dirWT, pattern=paste("chr"), full.names=TRUE)

system(paste0("rm comp_HiTC_sub_",as.character(binsize),"_",type,".bed"))
for(files in all_files){
  
  d=read.table(files,skip=1,row.names=1,header=T)
  for(i in 1:(as.integer(nrow(d)/nbin)+1)){
  
    write.table(x="#prova",file="prova",quote=F,col.names=F,row.names=F)
    write.table(x=d[((i-1)*nbin+1):min(i*nbin,nrow(d)),((i-1)*nbin+1):min(i*nbin,nrow(d))],file="prova",sep="\t",col.names=NA,quote=F,append=TRUE)
    system("sed -i '2,2s/\\./|/g' prova")
    lWT <- sapply("prova",
                  import.my5C)
    hiC_WT <- HTClist(lWT)

    pr1 <- pca.hic(hiC_WT[[1]], npc=2, asGRangesList=TRUE,gene.gr=H3K27ac_gr,logbin=TRUE)
    if(!is.null(pr1)){
      pr = pr1[[1]]
      
      scoreWT=c(mcols(pr)$score)
      colna=which(is.na(scoreWT))
      a=data.frame(chr=as.vector(seqnames(pr)),start=as.vector(start(pr)),end=as.vector(end(pr)),scoreWT=scoreWT)
      write.table(x=a,file=paste0("comp_HiTC_sub_",as.character(binsize),"_",type,".bed"),quote=F,col.names=F,row.names=F,append=T)
      hiC14norm <- normPerExpected(hiC_WT[[1]], method="mean")
      mapC(HTClist(hiC14norm), log.data=TRUE,tracks=GRangesList(H3K27ac_gr))
      }else{warning(paste(files," does not produce data"))}
  }
}

