rm(list=ls())
source("vars.R")
source("functions.R")

#################


### overall read coverage
if(T){
  load(paste0(DATADIR,"CallsForHeatMaps_gm12878_50xsc.RData"))
  gm50 <- mycalls(paste0(RAWDIR,"gm12878_50xsc.candidate.insertions.bed"),filt = T)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  
  ## melt call overlap
  m<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  m$V5 <- gsub("<INS:ME:","",m$V5)
  m$V5 <- gsub(">","",m$V5)
  m$V5 = gsub("ALU","Alu",m$V5)
  m$V5 = gsub("LINE1","L1",m$V5)
  m =  with(m,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  res$melt_pass <- countOverlaps(resize(res,101,"center"),m[m$Evi=="PASS"],ignore.strand=T)
  res$melt_pass <- ifelse(res$melt_pass>1,1,res$melt_pass)
  
  m<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  m[grep("0/1|1/1",m$V10),]$V7 <- "PASS"
  m$V5 <- gsub("<INS:ME:","",m$V5)
  m$V5 <- gsub(">","",m$V5)
  m$V5 = gsub("ALU","Alu",m$V5)
  m$V5 = gsub("LINE1","L1",m$V5)
  m =  with(m,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  res$melt_gt <- countOverlaps(resize(res,101,"center"),m[m$Evi=="PASS"],ignore.strand=T)
  res$melt_gt <- ifelse(res$melt_gt>1,1,res$melt_gt)
  
  d <- read.delim("C:/d/cov_comparison/GM12878_HiC_1200M.GenomeCov.txt",header=F)
  e <- read.delim("C:/d/cov_comparison/GM12878_WGS_1200M.GenomeCov.txt",header=F)
  d$id <- paste(d$V1,d$V2,d$V3,sep="_")
  e$id <- paste(e$V1,e$V2,e$V3,sep="_")
  names(d)[4] <- "HiC"
  names(e)[4] <- "WGS"
  e <- e[,c("WGS","id")]
  d <- d[,c("HiC","id")]
  d <- merge(d,e,by="id",all.x=T,all.y=T)
  d$chr <- gsub("_\\S*","",d$id)
  d$start <- gsub("chr\\S*?_","",d$id)
  d$start <- gsub("_\\S*","",d$start)
  d$end <- gsub(".*_","",d$id)
  d$start <- gsub("_\\S*","",d$start)
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  d <- d[!is.na(d$start),]
  rm(e)
  d <- with(d,GRanges(chr,IRanges(start,end),"+",HiC,WGS))
  d$HiC <- ifelse(is.na(d$HiC),0,d$HiC)
  d$WGS <- ifelse(is.na(d$WGS),0,d$WGS)
  
  d$HiC <- round(d$HiC*1e9/sum(d$HiC),2)
  d$WGS <- round(d$WGS*1e9/sum(d$WGS),2)
  d$class = cut(log2(d$WGS+0.1),breaks = 200)
  table(d$class)
  range(log2(d$WGS+0.1))
  
  o <- findOverlaps(res,d,ignore.strand=T)
  res$HiC <- res$WGS <- 0
  res[queryHits(o)]$HiC <- d[subjectHits(o)]$HiC
  res[queryHits(o)]$WGS <- d[subjectHits(o)]$WGS
  
  res <- resize(res,1)
  r <- res
  r$hitea <- 0
  r[r$call%in%c("Germline","HiTEA+PacB","unique")]$hitea <- 1
  table(r$call)
  r <- r[r$call%in%c("HiTEA+PacB","PacB Germline","PacB only")]
  r <- r[r$melt_gt>0 | r$hitea>0]
  
  r$l2f <- log2((r$HiC+.1)/(r$WGS+.1))
  z <- data.frame(class=r$melt_gt,hitea=r$hitea,WGS=r$WGS, HiC=r$HiC,change=r$l2f)
  z$class <- paste0(z$class,z$hitea)
  z$def=NA
  z[grep("01",z$class),]$def <- "158 missed by\n MELT (GT)"
  z[grep("10",z$class),]$def <- "183 missed by \nHiTEA"
  z[grep("11",z$class),]$def <- "948 correctly called by\nHiTEA and MELT (GT)"
  table(z$def)
  z <- z[z$class!="10",]
  save(z,res, file=paste0(DATADIR,"Coverage_CallsMissedbyMELTHiTEA_GM1287850x.RData"))
  
}

## where are pacb germline calls?
if(F){
  load("A:/work/hg38_annotations/Hg38_repeats.RData")
  load(paste0(DATADIR,"CallsForHeatMaps_gm12878_50xsc.RData"))
  gm50 <- mycalls(paste0(RAWDIR,"gm12878_50xsc.candidate.insertions.bed"),filt = T)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  
  ## Repeat annotation
  o = findOverlaps(gm50,hg38,ignore.strand=T)
  gm50$rep <- NA
  gm50[queryHits(o)]$rep <- as.character(hg38[subjectHits(o)]$repFamily)
  o = findOverlaps(res,hg38,ignore.strand=T)
  res$rep <- NA
  res[queryHits(o)]$rep <- as.character(hg38[subjectHits(o)]$repFamily)
  rm(o)
  
  ## melt call overlap
  m<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  m$V5 <- gsub("<INS:ME:","",m$V5)
  m$V5 <- gsub(">","",m$V5)
  m$V5 = gsub("ALU","Alu",m$V5)
  m$V5 = gsub("LINE1","L1",m$V5)
  m =  with(m,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  m <- m[m$TE%in%c("Alu","SVA","L1")]
  res$melt_pass <- countOverlaps(resize(res,101,"center"),m[m$Evi=="PASS"],ignore.strand=T)
  res$melt_pass <- ifelse(res$melt_pass>1,1,res$melt_pass)
  rm(m)
  m<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  m[grep("0/1|1/1",m$V10),]$V7 <- "PASS"
  m$V5 <- gsub("<INS:ME:","",m$V5)
  m$V5 <- gsub(">","",m$V5)
  m$V5 = gsub("ALU","Alu",m$V5)
  m$V5 = gsub("LINE1","L1",m$V5)
  m =  with(m,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  m <- m[m$TE%in%c("Alu","SVA","L1")]
  res$melt_gt <- countOverlaps(resize(res,101,"center"),m[m$Evi=="PASS"],ignore.strand=T)
  res$melt_gt <- ifelse(res$melt_gt>1,1,res$melt_gt)
  res$lowmelt <- countOverlaps(resize(res,101,"center"),m[m$Evi!="PASS"],ignore.strand=T)
  res$lowmelt <- ifelse(res$lowmelt>1,1,res$lowmelt)
  rm(m)
  
  
  ## calls overlapping to SVs in GM12878
  load("C:/d/gm12878_pacB/HG001_GRCh38_PBSV_Indels_usableSV.RData")
  #sv <- sv[sv$type%in%c("INS")]
  #sv <- sv[width(sv)>50]
  res$sv <- countOverlaps(resize(res,201,"center"),sv,ignore.strand=T)
  res$sv <- ifelse(res$sv>1,1,res$sv)
  o <- findOverlaps(resize(res,201,"center"),sv,ignore.strand=T)
  res$SVtype <- ""
  res[queryHits(o)]$SVtype <- sv[subjectHits(o)]$V6
  #res$SVSize <- "
  #res[queryHits(o)]$SVsize <- width(sv[subjectHits(o)])
  rm(o,sv)
  
  ## check if the discordent reads are enriched on both sides
  load(file="C:/d/WGS_MELTRAM.RData")
  wgs <- wgs[-grep("S|H",wgs$cigar)]
  p <- GRanges()
  for(f in c("Alu","L1","SVA")){
    cat(f,"\n")
    g <- wgs[wgs$te==f]
    r <- res[res$TE==f]
    r$up <- countOverlaps(flank(r,500),g[strand(g)=="+"],ignore.strand=T) 
    r$down <- countOverlaps(flank(r,500,start = F),g[strand(g)=="-"],ignore.strand=T) 
    p <- c(p,r)
    rm(g,r)
  }
  res <- p
  rm(p)
  res$RAMtest <- ifelse(res$down/(res$up+res$down)<0.01 | res$up/(res$up+res$down) <0.01|(res$up+res$down)<5,1,0)
  sum(res$RAMtest)
  table(res$RAMtest,res$melt_pass)
  rm(wgs)
  
  save(res,gm50,file=paste0(DATADIR,"CallsMissedAnnotations_GM1287850x.RData"))
  
}