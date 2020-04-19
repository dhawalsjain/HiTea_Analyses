
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(EnrichedHeatmap)
library(reshape2)
library(tidyr)
library(grid)
library(gridExtra)
library(data.table)
library(ggplot2)
library(ggrepel)
library(wesanderson)
library(viridis)
require(ggpubr)

###########################################################################
### functions
###########################################################################
mycalls <- function(path,filt=T){
  ins <- read.delim(path,header = F,comment.char = "#")
  names(ins) <- c("chr","start","end","id","score",'strand',"TE","EVI","Descr","remark") 
  if(filt==T){
    ins <- ins[ins$EVI>0,]
  }
  ins <- with(ins,GRanges(chr,IRanges(start,end),strand,id,TE,EVI,score,Descr,remark))
  return(ins)
}

mycalls.germlineAnnotations <- function(gr){
  ## based on trio
  p = rbind(read.delim("A:/work/RT in HiC/MELT/GM12878_trio_MELT/NA12891_hg38.txt",header=F,stringsAsFactors = F),
            read.delim("A:/work/RT in HiC/MELT/GM12878_trio_MELT/NA12892_hg38.txt",header=F,stringsAsFactors = F))
  p = with(p,GRanges(V1,IRanges(V2,V2),"*",TE=V5,Evi=V7,Descr=V8))
  p$TE <- gsub("<INS:ME:","",p$TE)
  p$TE <- gsub(">","",p$TE)
  p$TE = gsub("ALU","Alu",p$TE)
  p$TE = gsub("LINE1","L1",p$TE)
  
  d <- read.delim("A:/work/hg38_annotations/1000k_TEIns_hg38.vcf.txt",header=F,comment.char = "#")
  d$AF <- stringr::str_extract(d$V8,";AF=\\S*?;")
  d$AF <- gsub(";","",d$AF)
  d$AF <- gsub("AF=","",d$AF)
  d$AF <-as.numeric(d$AF)
  d$V5 <- gsub("<INS:ME:","",d$V5)
  d$V5 <- gsub(">","",d$V5)
  d <- with(d,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=AF))
  d$TE <- gsub("<INS:ME:","",d$TE)
  d$TE <- gsub(">","",d$TE)
  d$TE = gsub("ALU","Alu",d$TE)
  d$TE = gsub("LINE1","L1",d$TE)
  
  gr$GL <- countOverlaps(gr,resize(p,200,"center"),ignore.strand=T)
  gr$GLp <- countOverlaps(gr,resize(d[d$Evi>0.1],200,"center"),ignore.strand=T)
  gr$GLp <- gr$GL+gr$GLp
  
  gr$GL <- ifelse(gr$GL>1,1,gr$GL)
  gr$GLp <- ifelse(gr$GLp>1,1,gr$GLp)
  gr$GL = ifelse(gr$GL==0,"Somatic","Germline")
  gr$GLp = ifelse(gr$GLp==0,"Somatic","Germline")
  gr
}

mycalls.germlineAnnotations_gno <- function(gr){
  ## based on trio
  d1 <- read.delim("A:/work/hg38_annotations/gnomadV2_hg38_TE.bed",header=F)
  d1$V5 = gsub("ALU","Alu",d1$V5)
  d1$V5 = gsub("LINE1","L1",d1$V5)
  d1$AF <- as.numeric(gsub(";","",gsub("AF=","",stringr::str_extract(d1$V7,";AF=\\S*?;"))))
  d1$EAS_AF <- as.numeric(gsub(";","",gsub("EAS_AF=","",stringr::str_extract(d1$V7,"EAS_AF=\\S*?;"))))
  d1$EUR_AF <- as.numeric(gsub(";","",gsub("EUR_AF=","",stringr::str_extract(d1$V7,"EUR_AF=\\S*?;"))))
  d1$AFR_AF <- as.numeric(gsub(";","",gsub("AFR_AF=","",stringr::str_extract(d1$V7,"AFR_AF=\\S*?;"))))
  d1$AMR_AF <- as.numeric(gsub(";","",gsub("AMR_AF=","",stringr::str_extract(d1$V7,"AMR_AF=\\S*?;"))))
  d1$SAS_AF <- as.numeric(gsub(";","",gsub("SAS_AF=","",stringr::str_extract(d1$V7,"SAS_AF=\\S*?;"))))
  d1 <- with(d1,GRanges(V1,IRanges(V2,V2),"+",TE=V5,AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF))
  d1 <- d1[d1$AF>0.001]
  
  gr$GLg <- countOverlaps(gr,resize(d1,200,"center"),ignore.strand=T)
  
  gr$GLg <- ifelse(gr$GLg>1,1,gr$GLg)
  gr$GLg = ifelse(gr$GLg==0,"Somatic","Germline")
  gr
}

callSummary.pacbio <- function(path="C:/d/GM12878_HindIII_r1.candidate.insertions.bed",
                                  pacb,wgs=m,xtea,filt=NULL){
  ht <- mycalls(path)
  ht_unf <- mycalls(path,filt = F)
  ht <- ht_unf[ht_unf$EVI>0]
  bw=201
  wgs <- wgs[wgs$Evi=="PASS"]
  
  ht <- ht[ht$TE%in%c("Alu","L1","SVA")]
  ht_unf <- ht_unf[ht_unf$TE%in%c("Alu","L1","SVA")]
  wgs <- wgs[wgs$TE%in%c("Alu","L1","SVA")]
  
  if(!is.null(filt)){
    ht <- ht[ht$TE==filt]
    ht_unf <- ht_unf[ht_unf$TE==filt]
    pacb <- pacb[pacb$TE==filt]
  }
  pacb<- unique(pacb)
  pacb$call <- "PacB only"
  pacb$xtmlt <- countOverlaps(resize(pacb,bw,"center"), c(wgs[,0],xtea[,0]),ignore.strand=T)
  pacb$hitea <- countOverlaps(resize(pacb,bw,"center"),ht,ignore.strand=T)
  pacb$call <- ifelse(pacb$xtmlt>0,as.character("xTEA/MELT+PacB"),as.character("PacB only"))
  if(length(pacb[pacb$hitea>0]$call)>0){
    pacb[pacb$hitea>0]$call <- as.character("HiTEA+PacB")
  }
  table(pacb$call)
  
  ht$call <- "HiTEA only"
  ht$pacb <- countOverlaps(ht,resize(pacb,bw,"center"),ignore.strand=T)
  if(length(ht[ht$pacb>0]$call)>0){
    ht[ht$pacb>0]$call <- as.character("HiTEA+PacB")
  }
  table(ht$call)
  
  s <- ht[ht$call!="HiTEA only"]
  if(length(s)>0){
    s$call <- "HiTEA+PacB"
    s$xtmlt <- countOverlaps(s, resize(c(m[,0],xtea[,0]),bw,"center" ),ignore.strand=T)
    if(length(s[s$xtmlt>0]$call)>0){
      s[s$xtmlt>0]$call <- as.character("HiTEA+PacB+\nxTEA/MELT")
    }
  }
  table(s$call)
  
  a <- ht[ht$call=="HiTEA only"]
  a$call <- "unique"
  a$xtmlt <- countOverlaps(a, resize(c(m[,0],xtea[,0]),bw,"center" ),ignore.strand=T)
  a[a$xtmlt>0]$call <- as.character("HiTEA+xTEA/MELT")
  a <- mycalls.germlineAnnotations(a)
  if(length(a[a$GLp=="Germline" & a$call=="unique"]$call)>0){
    a[a$GLp=="Germline" & a$call=="unique"]$call <- "Germline"
  }
  table(a$call)
  table(a$GLp)
  
  
  q <- ht[ht$call=="HiTEA only"]
  q$call <- "unique"
  q <- mycalls.germlineAnnotations(q)
  if(length(q[q$GLp=="Germline"]$call)>0){
    q[q$GLp=="Germline"]$call <- "Germline"
  }
  table(q$call)
  
  res <- c(pacb[,c("TE","call")],q[,c("TE","call")])
  table(res$call)
  res$call <- ifelse(res$call=="xTEA/MELT+PacB",as.character("PacB only"),as.character(res$call))
  
  xx <- res[res$call=="PacB only"]
  xx <- mycalls.germlineAnnotations(xx)
  if(length(xx[xx$GLp=="Germline"]$call)>0){
    xx[xx$GLp=="Germline"]$call <- "PacB Germline"
  }
  table(xx$call)
  res <- c(res[res$call!="PacB only"],xx)
  table(res$call)
  return(res)
  
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xpacBHiT.RData")
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xFullpacBHiT.RData")
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xapacBBHiT.RData")
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xsbpacBBHiT.RData")
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_HindIII_r1pacBBHiT.RData")
  #save(res,file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_10xapacBBHiT.RData")
  
}

myspesen.function.pacbio <- function(pacb,path="A:/work/RT in HiC/PAPER/RAW/gm12878_50xsc.candidate.insertions.bed",gr=NULL,spsn=T){
  LEN=201
  TE="Alu"
  alu <- resize(pacb[pacb$TE==TE],LEN,"center")
  TE="L1"
  l1 <- resize(pacb[pacb$TE==TE],LEN,"center")
  TE="SVA"
  sva <- resize(pacb[pacb$TE==TE],LEN,"center")
  
  if(is.null(gr)){
    ht <- mycalls(path,filt = F)
    ht <- ht[ht$TE%in%c("Alu","SVA","L1")]
    #ht$EVI <- ifelse(ht$EVI==2,0,ht$EVI)
  }else{
    ht <- gr
  }
  pacb <- unique(pacb)
  
  ht <- resize(ht,LEN,"center")
  ht$alu <- countOverlaps(ht,alu,ignore.strand=T)
  ht$l1 <- countOverlaps(ht,l1,ignore.strand=T)
  ht$sva <- countOverlaps(ht,sva,ignore.strand=T)
  
  
  if(spsn == T){
    ####################### Specificity/ sensitivity
    # TP = true calls identified correctly
    # FP = false calls identified incorrectly
    # FN = true calls incorrectly identified
    # TN = false calls identified correctly (can't estimate here!)
    pl <- data.frame(
      TP=c(sum(ht[ht$EVI>0 & ht$TE=="Alu"]$alu>0),sum(ht[ht$EVI>0 & ht$TE=="L1"]$l1>0),sum(ht[ht$EVI>0 & ht$TE=="SVA"]$sva>0)),
      FP=c(sum(ht[ht$EVI>0 & ht$TE=="Alu"]$alu==0),sum(ht[ht$EVI>0 & ht$TE=="L1"]$l1==0),sum(ht[ht$EVI>0 & ht$TE=="SVA"]$sva==0)),
      TN=c(sum(ht[ht$EVI==0 & ht$TE=="Alu"]$alu==0),sum(ht[ht$EVI<1 & ht$TE=="L1"]$l1==0),sum(ht[ht$EVI<1 & ht$TE=="SVA"]$sva==0))
    )
    pl <- cbind(pl,FN=c((length(alu) - sum(ht[ht$EVI>0 & ht$TE=="Alu"]$alu>0)),
                        (length(l1) - sum(ht[ht$EVI>0 & ht$TE=="L1"]$l1>0)),
                        (length(sva) - sum(ht[ht$EVI>0 & ht$TE=="SVA"]$sva>0))
    )) 
    
    pl <- rbind(pl,apply(pl[1:3,],2,sum))
    pl <- cbind(TE=c("Alu","L1","SVA","Total"),pl)
    pl$TotalCalls <- length(ht[ht$EVI>0])
    pl$sensitivity <- round(pl$TP/(pl$TP+pl$FN),3)
    pl$specificity <- round(pl$TN/(pl$TN+pl$FP),3)
    pl$fdr <- round(pl$FP/(pl$FP+pl$TP),3)
    pl$precision <- round(pl$TP/(pl$FP+pl$TP),3)
    #pl$F1 <- round((2*pl$sensitivity*pl$precision)/(pl$sensitivity+pl$recall),3)
    pl$missrate <- round(pl$FN/(pl$FN+pl$TP),3)
    return(pl)
  }else{
    return(ht)
  }
  
}


