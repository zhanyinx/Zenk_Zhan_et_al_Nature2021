
library(HiTC)
library(plot3D)
library(RColorBrewer)
library(Rcpp)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggplot2)

sourceCpp("RCPP_functions/matrix_resize.cpp")
sourceCpp("RCPP_functions/vector_resize.cpp")


##global variables
type = "HP1"
dirHP1 <- "HP1/matrix/"  #directory containing Hi-C pro output matrix in the matrix format for HiTC 
dirWT <- "WT/matrix/"  #directory containing Hi-C pro output matrix in the matrix format for HiTC 
exclude = c("chr4","chrY")
size = 50 #rebinned matrix size
limit=5000000 ##limit of distance where we do not use to calculate compscore and saddle plots
binsize=10000 ##binsize of matrix
quant=0.999 ##normalised interactions above this quantile are set to the value of the quantile
H3K27ac_peaks <- "utils/H3K27ac_keep1.idr_0.01.filtered"
local_compartments <- "utils/comp_HiTC_sub_10000_Iovino_WT_all_july2019.bed"
##Telomers and centromers
tc_gr = import("utils/dm6_centromeres_telomeres",format="BED")
#####


###H3K27ac for A/B compartments
H3K27ac_peaks_data = read.table(H3K27ac_peaks)
H3K27ac_peaks_gr = makeGRangesFromDataFrame(H3K27ac_peaks_data,
                                            seqnames.field="V1",
                                            start.field="V2",
                                            end.field="V3",
                                            strand.field="V4")
#End H3K27ac



##reading matrix file for HiTC
files=list.files(dirWT, pattern=paste("chr"), full.names=TRUE)
for(i in 1:length(exclude)){
  if(length(grep(exclude[i],files))>0){
    files = files[-grep(exclude[i],files)]
  }
}
lWT <- sapply(files,
              import.my5C)
hiC_WT <- HTClist(lWT)


files = list.files(dirHP1, pattern=paste("chr"), full.names=TRUE)
for(i in 1:length(exclude)){
  if(length(grep(exclude[i],files))>0){
    files = files[-grep(exclude[i],files)]
  }
}
lHP1 <- sapply(files,
                import.my5C)
hiC_HP1 <- HTClist(lHP1)

##Reading compartments
lcomp <- read.table(local_compartments)

##Initialiting variables
eigenvectorWT=NULL
eigenvectorHP1=NULL
allscores=NULL

## Performed PCA
contaWT=0
contaHP1_WTComp=0
contaHP1=0




for(i in 1:length(hiC_WT)){
  pr1<-pca.hic(hiC_WT[[i]], npc=2, asGRangesList=TRUE,gene.gr=H3K27ac_peaks_gr,logbin=TRUE) ##PCA
  pr = pr1[[1]]
  compartmentsWT = (merge(lcomp,data.frame(pr),by=c(1,2,3),sort=FALSE,all.y=TRUE))$V4 ##
  
  pr1<-pca.hic(hiC_HP1[[i]], npc=2, asGRangesList=TRUE,gene.gr=H3K27ac_peaks_gr,logbin=TRUE) ##PCA
  prHP1 = pr1[[1]]
  
  ###Extracting scores
  scoreWT=c(mcols(pr)$score) ##eigenvector
  colnaWT=which(is.na(scoreWT)) ##which are na
  scoreHP1=c(mcols(prHP1)$score) ##eigenvector
  colnaHP1=which(is.na(scoreHP1)) ##which are na
  
  ##removing telomers and centromers
  ov <- findOverlaps(subject=pr,query=tc_gr)
  
  colna=unique(c(colnaWT,colnaHP1,subjectHits(ov)))
  
  scoreWT=scoreWT[-colna] ##removing NA
  compartmentsWT=compartmentsWT[-colna] ##removing NA
  
  
  ##Extracting normalised data (distance normalised)
  norm = intdata(normPerExpected(hiC_WT[[i]], method="mean",logbin=TRUE)) 
  norm = norm[-colna,-colna] ##Removing NA row and columns
  norm[norm>quantile(as.matrix(norm),0.999)]=quantile(as.matrix(norm),0.999)
  
  ##limit on distance to calculate saddle plot and compscore
  appo=melt(as.matrix(norm))
  appo[,1]=as.numeric(appo[,1])
  appo[,2]=as.numeric(appo[,2])
  appo[abs(appo[,1]-appo[,2])*binsize>limit,3]=-1 ##IMPORTANT I PUT -1 to exclude only too long range interactions
  norm=acast(appo,Var1~Var2)
  
  
  
  ##Ordering matrix based on increasing eigenvector
  norm_appo=norm
  norm = norm[order(compartmentsWT),order(compartmentsWT)] 
  colnames(norm)=as.character(1:ncol(norm))
  rownames(norm)=as.character(1:ncol(norm))
  
  #Resizing matrix
  if(nrow(norm)>=size){
    appo = cresize_nozeros(as.matrix(norm),size)
   
    if(contaWT==0){
      binnedWT=appo
      contaWT=contaWT+1
    }else{
      binnedWT=binnedWT+appo
      contaWT=contaWT+1
    }  
  }  
  
  
  #HP1 using WT eigenvector
  norm = intdata(normPerExpected(hiC_HP1[[i]], method="mean",logbin=TRUE))
  norm = norm[-colna,-colna]
  if(length(compartmentsWT)!=nrow(norm)){
    warning("dimension of a matrix in HP1 is different from WT! Removing the exceeding parts")
  }
  norm[is.na(norm)]=0
  norm[norm>quantile(as.matrix(norm),0.999)]=quantile(as.matrix(norm),0.999)
  
  
  ##limit on distance to calculate saddle plot and compscore
  appo=melt(as.matrix(norm))
  appo[,1]=as.numeric(appo[,1])
  appo[,2]=as.numeric(appo[,2])
  appo[abs(appo[,1]-appo[,2])*binsize>limit,3]=-1
  norm=acast(appo,Var1~Var2)
  
 
  norm_appo=norm
  norm = norm[order(compartmentsWT),order(compartmentsWT)]
  if(nrow(norm)>=size){
    appo = cresize_nozeros(as.matrix(norm),size)
    if(contaHP1_WTComp==0){
      binnedHP1_WTcomp=appo
      contaHP1_WTComp=contaHP1_WTComp+1
    }else{
      binnedHP1_WTcomp=binnedHP1_WTcomp+appo
      contaHP1_WTComp=contaHP1_WTComp+1
    }  
  }
  
  
 
}


##Normalising matrix 
binnedWT=binnedWT/contaWT
binnedHP1_WTcomp=binnedHP1_WTcomp/contaHP1_WTComp






