rm(list=ls())
source("functions.R")
source("vars.R")

###########################################################################
## test 
if(F){
  myf<- function(pacb){
    h <- distanceToNearest(pacb) %>% as.data.frame()
    h$ref <- pacb[h$queryHits]$TE
    h$test <- pacb[h$subjectHits]$TE
    h <- h[h$distance<100,]
    h <- h[h$ref==h$test,]
    h$queryHits %>% unique
    h <- h[duplicated(t(apply(h[,1:2],1,sort))),]
    if(nrow(h)>0){
      pacb <- pacb[-h$queryHits]
    }
    pacb
  }
  getGR <- function(pacb, path="C:/d/gm12878sb.candidate.insertions.bed"){
    a <- myspesen.function.pacbio(pacb=pacb,path =path,spsn = F)
    a <- mycalls.germlineAnnotations(a)
    a$pacb <- a$alu+a$l1+a$sva
    a$pacb <- ifelse(a$pacb>1,1,a$pacb)
    a <- resize(a,1,"center")
    a$x <- paste0(seqnames(a),"_",start(a))
    a
  }
  pacb = rbind(read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_LINE1.txt",header=F),
               read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_Alu.txt",header=F),
               read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_SVA.txt",header=F))
  pacb <- with(pacb,GRanges(V1,IRanges(V2,V2),"*",TE=V3))
  pacb$TE <- gsub("LINE1","L1",pacb$TE)
  alu <- myf(unique(pacb[pacb$TE=="Alu"]))
  sva <- myf(unique(pacb[pacb$TE=="SVA"]))
  l1 <- myf(unique(pacb[pacb$TE=="L1"]))
  pacb<- c(alu,l1,sva)
  table(pacb$TE)
  rm(alu,sva,l1)
  
  myspesen.function.pacbio(pacb=pacb,path ="C:/d/gm12878sb.candidate.insertions.bed",spsn = T)
  myspesen.function.pacbio(pacb=pacb,path ="C:/d/gm12878_50xsb.candidate.insertions.bed",spsn = T)
  myspesen.function.pacbio(pacb=pacb,path ="C:/d/gm12878_20xsb.candidate.insertions.bed",spsn = T)
  
  
  myspesen.function.pacbio(pacb=pacb,path ="A:/work/RT in HiC/PAPER/RAW/gm12878_10xsc.candidate.insertions.bed",spsn = T)
  myspesen.function.pacbio(pacb=pacb,path ="A:/work/RT in HiC/PAPER/RAW/gm12878_20xsc.candidate.insertions.bed",spsn = T)
  myspesen.function.pacbio(pacb=pacb,path ="A:/work/RT in HiC/PAPER/RAW/gm12878_50xsc.candidate.insertions.bed",spsn = T)
  myspesen.function.pacbio(pacb=pacb,path ="A:/work/RT in HiC/PAPER/RAW/gm12878sc.candidate.insertions.bed",spsn = T)
  
  
  a <- getGR(pacb,"C:/d/gm12878sb.candidate.insertions.bed") 
  ax <- getGR(pacb,"A:/work/RT in HiC/PAPER/RAW/gm12878sc.candidate.insertions.bed")
  
  table(ax$TE,ax$l1)
  #ax$RAM <- gsub(";.*","",gsub(".*RAM=","",ax$Descr)) %>% as.numeric
  #ax$test <- round(ax$RAM/ax$score,2)
  ax$rec=gsub(";.*","",gsub(".*RECI%=","",ax$Descr))
  ax$ori=gsub(";.*","",gsub(".*isPolyA=","",ax$Descr))
  
  ### missed calls
  gx <- ax[ax$TE=="Alu" & ax$EVI==0 & ax$alu>0] %>% as.data.frame
  gx <- ax[ax$TE=="L1" & ax$EVI==0 & ax$l1>0] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI==0 & ax$sva>0] %>% as.data.frame
  
  gx1 <- gx[,c(1:3,6,7,5)]
  gx1$strand <- "+"
  write.table(gx1,file="C:/d/vis/checks.bed",quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## FP
  gx <- ax[ax$TE=="L1" & ax$EVI>0 & ax$l1==0] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI>0 & ax$sva==0] %>% as.data.frame
  gx <- ax[ax$TE=="Alu" & ax$EVI>0 & ax$alu==0 & ax$GLp=="Somatic"] %>% as.data.frame
  
  ## TP
  gx <- ax[ax$TE=="L1" & ax$EVI>0 & ax$l1>0] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI>0 & ax$sva>0] %>% as.data.frame
  
  ### chk: L1
  ### 174315623; 21130729
  
  
  
  g <- a[a$EVI>0 & a$pacb==0 & a$GLp=="Somatic"] %>% as.data.frame
  gx <- ax[ax$EVI>0 & ax$pacb==0 & ax$GLp=="Somatic"] %>% as.data.frame
  gx <- gx[!gx$x%in%g$x,]
  
  g <- a[a$TE=="Alu" & a$EVI>0 & a$alu==0 & a$GL=="Somatic"] %>% as.data.frame
  gx <- ax[ax$TE=="Alu" & ax$EVI>0 & ax$alu==0 & ax$GL=="Somatic"] %>% as.data.frame
  
  g <- a[a$TE=="SVA" & a$EVI>0 & a$sva==0 ] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI>0 & ax$sva==0 ] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI==0 & ax$sva>0 ] %>% as.data.frame
  
  g <- a[a$TE=="L1" & a$EVI>0 & a$l1==0 & a$GL=="Somatic"] %>% as.data.frame
  gx <- ax[ax$TE=="L1" & ax$EVI>0 & ax$l1==0 & ax$GL=="Somatic"] %>% as.data.frame
  
  g <- a[a$TE=="L1" & a$EVI>0 & a$l1>0 ] %>% as.data.frame
  gx <- ax[ax$TE=="L1" & ax$EVI>0 & ax$l1>0] %>% as.data.frame
  g <- g[!g$x%in%gx$x,]
  
  g <- a[a$TE=="L1" & a$EVI>0 & a$l1==0 ] %>% as.data.frame
  gx <- ax[ax$TE=="L1" & ax$EVI>0 & ax$l1==0] %>% as.data.frame
  gx <- gx[!gx$x%in%g$x,]
  a[a$x%in%gx$x]
  
  g <- a[a$TE=="SVA" & a$EVI>0 & a$sva==0 ] %>% as.data.frame
  gx <- ax[ax$TE=="SVA" & ax$EVI>0 & ax$sva==0] %>% as.data.frame
  gx <- gx[!gx$x%in%g$x,]
  a[a$x%in%gx$x]
  

} 


### beta binomial for disc reads
if(F){
  disc <- read.delim("C:/d/temp/20x.disc.txt.gz",header=F)
  disc <- with(disc,GRanges(V1,IRanges(V2,V2),"*",all=V3,k1=V4,k2=V5,k5=V6,k10=V7))
  
  a <- mycalls(path = "C:/d/temp/gm12878_20xsb.candidate.insertions.bed",filt = F)
  a$name <- a$id
  p <- read.delim("C:/d/temp/gm12878_20xsb_RandomLocs.bed",header=F,comment.char = "#")
  names(p) <- c("chr","start","end","name","score","strand","cov")
  p$rname <- paste0("random_",1:length(p$chr))
  p <- as(p,"GRanges")
  p <- resize(p,1)
  p$cov <- p$cov*10
  
  load("A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xsbpacBBHiT.RData")
  res$name <- 1:length(res)
  
  myf <- function(p,disc){
    o <- findOverlaps(resize(p,501,"center"),disc,ignore.strand=T) %>% as.data.frame
    o$id <- p[o$queryHits]$name
    o <- cbind(o,as.data.frame(disc[o$subjectHits]))
    o <- o[,c("id","all","k1","k2","k5","k10")]
    o <- data.table(o)
    o <- o[, lapply(.SD, sum, na.rm=TRUE), by=id ]
    o <- as.data.frame(o)
    o$rk1 <-  (o$all-o$k1)/o$all
    o$rk2 <-  (o$all-o$k2)/o$all
    o$rk5 <-  (o$all-o$k5)/o$all
    o$rk10 <-  (o$all-o$k10)/o$all
    return(o)
  }
  
  pf <- myf(p,disc)
  af <- myf(a,disc)
  rf <- myf(res,disc)
  
  hist(pf$rk10,breaks=100,col="gray",probability = T)
  hist(af$rk10,breaks=100,add=T,col="red",probability = T)
  hist(rf$rk10,breaks=100,add=T,col="green",probability = T)
  
  median(pf$rk10)
  median(af$rk10)
  median(rf$rk10)
  
  par(mfrow=c(2,2))
  smoothScatter(log(pf$all),pf$rk10)
  smoothScatter(log(af$all),af$rk10)
  smoothScatter(log(rf$all),rf$rk10)
  
  
  
}




## plot tests
if(F){
  ## 1
  nrows <- length(plotdf)
  
  par(mfrow=c(nrows,3))
  for( f in names(plotdf)){
    pl <- plotdf[[paste0(f)]]
    
    ## smooth scatter plot
    smoothScatter( pl[pl$dataframe=="ctr",]$x, pl[pl$dataframe=="ctr",]$y,xlab="Total coverage",ylab="Hi-C ambiguous read coverage",cex.lab=1.5,cex.axis=1.5,cex.main=3)
    loess.fit = tryCatch({
      loess.smooth(pl[pl$dataframe=="ctr",]$x, pl[pl$dataframe=="ctr",]$y)
    },error = function(e) { return(NA)})
    if(!is.na(loess.fit)){
      lines(loess.fit,col="red")
    }
    ## mead vs sd plot    
    plot(pl[pl$dataframe=="bg",]$x,pl[pl$dataframe=="bg",]$y,main=paste0(f),xlab="mean",ylab="sd",cex=1.2,cex.lab=1.5,cex.axis=1.5,cex.main=3,pch=16,col=rgb(1,0,0,0.2))
    lines(x = pl[pl$dataframe=="bg",]$x,y=sqrt(pl[pl$dataframe=="bg",]$x),col="black")
    
    ## P-value histogram  
    hist(pl[pl$dataframe=="ex",]$x,breaks=20,xlab="p-value",ylab = "frequency",main=paste0(f," (p-values)"),cex.lab=1.5,cex.axis=1.5,cex.main=3)
  }
  
  
  
  ## 2
  file=paste0(params$dir,"/",params$outprefix,".ReadPairInsertSizeOriSummary.logs.gz")
  pl <- read.delim(file = paste0(file),header = F)
  
  ggplot(pl,aes(x=V1,y=log10(V3),col=V2))+geom_line(size=1.5)+
    theme_bw()+ xlab("inter-mate distance, log10")+ ylab("number of reads, log10")+
    geom_vline(xintercept = c(3,3.3),col="gray",linetype="dashed")+
    theme(axis.title = element_text(size=18,color="black"),
          axis.text = element_text(size=15,colour = "black"),
          legend.position = "top",
          legend.text = element_text(size=16,colour = "black"),
          legend.title = element_blank())
  
  
  
  #x <- paste0(rep("N",171823),collapse = "")
  #writeXStringSet(DNAStringSet(x,start = 1),file=paste0(DATADIR,"chrEBV.fa"))
  #d <- 
  #rbind(  
  #cbind(read.delim("A:/work/scripts/TE_insertions/final/v_paper/hg19/ALU.bed",header = F,stringsAsFactors = F),TE="Alu"),
  #cbind(read.delim("A:/work/scripts/TE_insertions/final/v_paper/hg19/LINE1.bed",header = F,stringsAsFactors = F),TE="L1"),
  #cbind(read.delim("A:/work/scripts/TE_insertions/final/v_paper/hg19/SVA.bed",header = F,stringsAsFactors = F),TE="SVA")
  #)
  #d$V1 <- paste0("chr",d$V1)
  #write.table(d,file="A:/work/scripts/TE_insertions/final/v_paper/hg19/bgRepeats_hg19.bed",quote = F,sep = "\t",row.names = F,col.names = F)

  d <- read.delim("A:/work/scripts/TE_insertions/final/v_paper/hg19/bgRepeats_hg19.bed",header = F,stringsAsFactors = F)
  d <- d[d$V7!="HERVK",]
  d$V5 <- "."
  #d$V4 <- "."
  write.table(d,file="A:/work/scripts/TE_insertions/final/v_paper/hg19/bgRepeats_hg19.bed",quote = F,sep = "\t",row.names = F,col.names = F)
  
  
    
  d <- read.delim("A:/work/scripts/TE_insertions/final/v_paper/hg38/bgRepeats_hg38.bed",header = F,stringsAsFactors = F)
  d <- d[d$V7!="HERVK",]
  d$V5 <- "."
  #d$V4 <- "."
  write.table(d,file="A:/work/scripts/TE_insertions/final/v_paper/hg38/bgRepeats_hg38.bed",quote = F,sep = "\t",row.names = F,col.names = F)
  
}