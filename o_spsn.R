rm(list=ls())
source("functions.R")
source("vars.R")

########################################################################################3
### get calls
if(T){
  
  ### xtea by simon
  xtea <- rbind(cbind(read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_illumina_Alu.txt",header=F),TE="Alu"),
                cbind(read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_illumina_L1.txt",header=F),TE="L1"),
                cbind(read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_illumina_SVA.txt",header=F),TE="SVA"))
  xtea <- unique(xtea[,c(1,2,50)])
  xtea <- with(xtea,GRanges(V1,IRanges(V2,V2),"+",TE=TE))
  
  ## MELT-GM12878
  m<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  m$V5 <- gsub("<INS:ME:","",m$V5)
  m$V5 <- gsub(">","",m$V5)
  m$V5 = gsub("ALU","Alu",m$V5)
  m$V5 = gsub("LINE1","L1",m$V5)
  m =  with(m,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  
  ## MELT-GM12878 (at 600M reads)
  m0<- read.delim("A:/work/RT in HiC/melt/gm12878_600M_5m/Insertions.txt",comment.char = "#",header=F)
  m0$V5 <- gsub("<INS:ME:","",m0$V5)
  m0$V5 <- gsub(">","",m0$V5)
  m0$V5 = gsub("ALU","Alu",m0$V5)
  m0$V5 = gsub("LINE1","L1",m0$V5)
  m0 =  with(m0,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  
  ## MELT-K562
  n<- read.delim("A:/work/RT in HiC/MELT/k562_5m/Insertions.txt",comment.char = "#",header=F)
  n$V5 <- gsub("<INS:ME:","",n$V5)
  n$V5 <- gsub(">","",n$V5)
  n$V5 = gsub("ALU","Alu",n$V5)
  n$V5 = gsub("LINE1","L1",n$V5)
  n =  with(n,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  
  ## pacb by simon
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
  
}

if(F){
  pacb <- as.data.frame(pacb)
  rd <- read.delim("C:/d/bgmodel/gm12878_10xsb_RandomLocs.bed",header = F)
  rd$TE = "random"
  names(rd) <- c("seqnames","start","end","width","d","strand","d1","TE")
  rd$d <- rd$d1 <- NULL
  rd <- rd[1:1000,]  
  pacb <- rbind(pacb,rd)
  write.table(pacb,file = paste0(DATADIR,"pacb_references_hg38.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
}


### 
if(F){
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
  
  d <- data.table(data.frame(d))
  d$Evi <- round(d$Evi,2)
  d$test <- 1
  d <- d[,tot:=sum(test),by=list(Evi,TE)]
  d <- as.data.frame(d)
  d1 <- d[order(d$Evi,decreasing = T),]
  d1 <- d1[,c("TE","Evi","tot")] %>% unique
  
  rect1=data.frame(xmin=0.1,xmax=Inf,ymin=0,ymax=Inf)
  
  
  ggplot(d1,aes(Evi,tot,col=TE))+ geom_point(size=3)+stat_smooth()+ scale_y_log10()+#stat_ecdf(geom = "point",size=3)+
    xlab("Population allele fraction")+ylab("number of total insertion calls")+
    geom_rect(data=rect1, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="gray", alpha=0.1,inherit.aes = F)+
    scale_x_reverse()
    #scale_x_continuous(expand = c(0,0))#+scale_y_continuous(expand = c(0,0))
  
  
  
}


### calls for display in the paper
if(F){
  for(f in c("gm12878_50xsc","gm12878_20xsc",'gm12878_10xsc','GM12878_HindIII_r1c')){
    cat(f,"\n")
    res <- callSummary.pacbio(path=paste0(RAWDIR,f,".candidate.insertions.bed"),
                              pacb=pacb,wgs=m,xtea,filt=NULL)
    save(res,file=paste0(DATADIR,"CallsForHeatMaps_",f,".RData"))
    res <- resize(res,width(res)+10001,"center")
    res <- as.data.frame(res)
    res$id <- 1:nrow(res)
    res <- res[,c(1:3,10,4,5)]
    write.table(res,file=paste0(DATADIR,f,"CallForHeatMaps.bed"),sep = "\t",quote = F,row.names = F,col.names = F)
  }
  
  load(file=paste0(DATADIR,"CallsForHeatMaps_gm12878_50xsc.RData"))
  load(file=paste0(DATADIR,"CallsForHeatMaps_GM12878_HindIII_r1c.RData"))
  
  table(res$call)
  res <- as.data.frame(res)
}



### TP/FP rates
if(T){
  res <- c()
  for(f in c("gm12878_10xsc","gm12878_20xsc","gm12878_50xsc","gm12878sc",
             "GM12878_HindIII_r1c","GM12878_HindIII_r2c","GM12878_HindIII_r3c","gm12878_ncoIc")){
    path=paste0(RAWDIR,f,".candidate.insertions.bed")
    cat(f,"\n",path,"\n\n")
    res <- rbind(res,cbind(myspesen.function.pacbio(pacb,path),profile=f))
  }
  
  gr <- m
  gr$EVI <- ifelse(gr$Evi==as.character("PASS"),3,0)
  gr <- gr[gr$TE%in%c("Alu","SVA","L1")]
  res <- rbind(res,cbind(myspesen.function.pacbio(pacb,"",gr),profile="WGS_MELT"))
  
  n<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F,stringsAsFactors = F)
  n[grep("0/1|1/1|1/0",n$V10),]$V7 <- "PASS"
  n$V5 <- gsub("<INS:ME:","",n$V5)
  n$V5 <- gsub(">","",n$V5)
  n$V5 = gsub("ALU","Alu",n$V5)
  n$V5 = gsub("LINE1","L1",n$V5)
  n =  with(n,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  gr <- n
  gr$EVI <- ifelse(gr$Evi==as.character("PASS"),3,0)
  gr <- gr[gr$TE%in%c("Alu","SVA","L1")]
  res <- rbind(res,cbind(myspesen.function.pacbio(pacb,"",gr),profile="WGS_MELT_Alt"))
  
  gr <- m0
  gr$EVI <- ifelse(gr$Evi==as.character("PASS"),3,0)
  gr <- gr[gr$TE%in%c("Alu","SVA","L1")]
  res <- rbind(res,cbind(myspesen.function.pacbio(pacb,"",gr),profile="WGS_MELT_600M"))
  
  n<- read.delim("A:/work/RT in HiC/melt/gm12878_600M_5m/Insertions.txt",comment.char = "#",header=F,stringsAsFactors = F)
  n[grep("0/1|1/1|1/0",n$V10),]$V7 <- "PASS"
  n$V5 <- gsub("<INS:ME:","",n$V5)
  n$V5 <- gsub(">","",n$V5)
  n$V5 = gsub("ALU","Alu",n$V5)
  n$V5 = gsub("LINE1","L1",n$V5)
  n =  with(n,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  gr <- n
  gr$EVI <- ifelse(gr$Evi==as.character("PASS"),3,0)
  gr <- gr[gr$TE%in%c("Alu","SVA","L1")]
  res <- rbind(res,cbind(myspesen.function.pacbio(pacb,"",gr),profile="WGS_MELT_600M_Alt"))
  
  
  gr <- xtea
  gr$EVI=3
  gr <- gr[gr$TE%in%c("Alu","SVA","L1")]
  #res <- rbind(res,cbind(myspesen.function(m,xtea,l1,Alu,"",gr),profile="WGS_xTEA"))
  res <- rbind(res,cbind(myspesen.function.pacbio(pacb,"",gr),profile="WGS_xTEA"))
  
  
  res$profile <- gsub("gm12878_10xsc","300M (MboI)",res$profile)
  res$profile <- gsub("gm12878_20xsc","600M (MboI)",res$profile)
  res$profile <- gsub("gm12878_50xsc","1.4B (MboI)",res$profile)
  res$profile <- gsub("gm12878sc","5B (MboI)",res$profile)
  res$profile <- gsub("GM12878_HindIII_r1c","1.8B (HindIII)",res$profile)
  res$profile <- gsub("GM12878_HindIII_r2c","670M (HindIII)",res$profile)
  res$profile <- gsub("GM12878_HindIII_r3c","1.9B (HindIII)",res$profile)
  res$profile <- gsub("gm12878_ncoIc","670M (NcoI)",res$profile)
  res$TE <- gsub("Total","All insertions",res$TE)
  res$F1 <- round(2*res$sensitivity*res$precision/(res$sensitivity+res$precision),3)
  res$recall <- res$sensitivity
  save(res,file=paste0(DATADIR,"SPSN_dataframe.RData"))
  
}
### 