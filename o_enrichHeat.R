rm(list=ls())
source("functions.R")
source("vars.R")

#######
if(F){
    hi <- read.delim(paste0(RAWDIR,"gm12878_50xsc.grfile.txt"),header=F)
    names(hi) <- c("chr","start","end","strand","evi","TE","mm")
    hi <- with(hi,GRanges(chr,IRanges(start,end),strand,TE,evi,mm))
    hi <- resize(hi,1,"start")
    save(hi,file=paste0(DATADIR,"gm50xsc_readSummary.RData"))
    
    library(ShortRead)
    wgs <- readGAlignments(paste0(RAWDIR,"WGS_50x.vis.disc.bam"))
    wgs <- with(wgs,GRanges(seqnames,IRanges(start,end),strand,cigar))
    wgs <- resize(wgs,1,"start")
    save(wgs,file=paste0(DATADIR,"wgs_discreadSummary.RData"))
}
  
################################## 50X data
if(F){
  load("C:/d/enrichheatmaps/wgs_readSummary.RData")
  load("C:/d/enrichheatmaps/gm12878_50xsb_readSummary.RData")
  hi$TE <- gsub(":\\S*","",hi$TE)
  hide$TE <- gsub(":\\S*","",hide$TE)
  wgsde <- wgs[grep("DE,TP",wgs$evi)]
  wgs <- wgs[-grep("DE,TP",wgs$evi)]
  #hide <- hi[grep("DE,TP",hi$evi)]
  #hi <- hi[-grep("DE,TP",hi$evi)]
  if(F){ ## proportions of Alu, L1, SVAs in the data
    hide$TE <- gsub("~\\S*","",hide$TE)
    hide$TE <- gsub("AluYa5","Alu",hide$TE)
    table(hide$TE)/length(hide)
    
    hi$TE <- gsub("~\\S*","",hi$TE)
    hi$TE <- gsub("AluYa5","Alu",hi$TE)
    table(hi$TE)/length(hi)
    
    
    table(wgs$TE)/length(wgs)
    table(wgsde$TE)/length(wgsde)
    
  }
  if(F){
    hi$TE <- hide$TE <- wgsde$TE <- wgs$TE <- NULL
    hi <- unique(hi)
    hide <- unique(hide)
    wgsde <- unique(wgsde)
    wgs <- unique(wgs)
    
  }
  
  
  load(paste0(DATADIR,"CallsForHeatMaps_gm12878_50xsc.RData"))
  table(res$call)
  res$id <- paste0("id",1:length(res))
  m.spl.wgs <- c()
  m.unspl.wgs <- c()
  m.unspl.hi <- c()
  m.spl.hi <- c()
  
  TE="Alu"
  win=4010
  bw=10
  for(TE in c("Alu","SVA","L1")){
    cat(TE,"\n")
    a1 <-res[res$TE==TE]
    ind <- resize(a1,win,"center")
    ind <- unlist(tile(ind,width = bw))
    
    gr1 <- wgs[wgs$TE%in%c(TE,"PolyA")] #,"PolyA"
    gr2 <- hi[hi$TE%in%c(TE,"PolyA")] #,"PolyA"
    #if(TE=="L1"){
    #  gr1 <- gr1[gr1$mm<3]
    #  gr2 <- gr2[gr2$mm<3]
    #}
    
    ind$wgs <- countOverlaps(ind,gr1,ignore.strand=T)
    ind$wgsde <- countOverlaps(ind,wgsde[wgsde$TE==TE],ignore.strand=T) # | wgsde$TE=="PolyA"
    ind$hi <- countOverlaps(ind,gr2,ignore.strand=T)
    ind$hide <- countOverlaps(ind,hide[hide$TE==TE],ignore.strand=T) # | hide$TE=="PolyA"
    
    x <- as.data.frame(matrix(ind$wgsde,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    y <- as.data.frame(matrix(ind$wgs,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    x1 <- as.data.frame(matrix(ind$hide,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    y1 <- as.data.frame(matrix(ind$hi,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    rownames(x) <- rownames(y) <- rownames(x1) <- rownames(y1) <- a1$id
    
    m.unspl.wgs <- rbind(m.unspl.wgs,y)
    m.spl.wgs <- rbind(m.spl.wgs,x)
    m.unspl.hi <- rbind(m.unspl.hi,y1)
    m.spl.hi <- rbind(m.spl.hi,x1)
    
    rm(x,y,x1,y1,ind,a1,gr1,gr2)
  }
  rm(wgs,hi,wgsde,hide)
  m.spl.hi <- m.spl.hi[match(res$id,rownames(m.spl.hi)),]
  m.unspl.hi <- m.unspl.hi[match(res$id,rownames(m.unspl.hi)),]
  m.spl.wgs <- m.spl.wgs[match(res$id,rownames(m.spl.wgs)),]
  m.unspl.wgs <- m.unspl.wgs[match(res$id,rownames(m.unspl.wgs)),]
  save(res,m.spl.hi,m.spl.wgs,m.unspl.hi,m.unspl.wgs,file=paste0(DATADIR,"CoveragePlotMatrices_50Xsc.RData"))
  
}


#### plotting WGS dis reads as called by MELT
if(F){
  load("C:/d/enrichheatmaps/wgs_discreadSummary.RData")
  load(paste0(DATADIR,"CallsForHeatMaps_gm12878_50xsc.RData"))
  table(res$call)
  res$id <- paste0("id",1:length(res))
  wgsde <- wgs[grep("S|H",wgs$cigar)]
  wgs <- wgs[-grep("S|H",wgs$cigar)]
  
  m.spl.wgs <- c()
  m.unspl.wgs <- c()
  win=4010
  bw=10
  a1 <-res
  ind <- resize(a1,win,"center")
  ind <- unlist(tile(ind,width = bw))
  ind$wgs <- countOverlaps(ind,wgs,ignore.strand=T)
  ind$wgsde <- countOverlaps(ind,wgsde,ignore.strand=T) # | wgsde$TE=="PolyA"
  m.spl.wgs <- as.data.frame(matrix(ind$wgsde,nrow = length(a1),ncol = round(win/bw), byrow = T) )
  m.unspl.wgs <- as.data.frame(matrix(ind$wgs,nrow = length(a1),ncol = round(win/bw), byrow = T) )
  rownames(m.spl.wgs) <- rownames(m.unspl.wgs) <- a1$id
  save(res,m.spl.wgs,m.unspl.wgs,file=paste0(DATADIR,"CoveragePlotMatrices_50Xsc_discWGS.RData"))
  
}


## numbers
if(F){ ## proportion of WGS-MELT reads
  library(ShortRead)
  alu <- readGAlignments("C:/d/vis/GM12878_hg38sorted.ALU.pulled.sorted.bam")
  l1 <- readGAlignments("C:/d/vis/GM12878_hg38sorted.LINE1.pulled.sorted.bam")
  sva <- readGAlignments("C:/d/vis/GM12878_hg38sorted.SVA.pulled.sorted.bam")
  
  alu <- with(alu,GRanges(seqnames(alu),IRanges(start(alu),end(alu)),strand(alu),cigar))
  alu$TE <- "Alu"
  l1 <- with(l1,GRanges(seqnames(l1),IRanges(start(l1),end(l1)),strand(l1),cigar))
  l1$TE <- "L1"
  sva <- with(sva,GRanges(seqnames(sva),IRanges(start(sva),end(sva)),strand(sva),cigar))
  sva$TE <- "SVA"
  wgs <- c(alu,l1,sva)
  rm(alu,l1,sva)
  wgsde <- wgs[grep("S|H",wgs$cigar)]
  wgs <- wgs[-grep("S|H",wgs$cigar)]
  
  table(wgsde$TE)/length(wgsde)
  table(wgs$TE)/length(wgs)
  
}


############## 20X data
if(F){
  load("C:/d/enrichheatmaps/gm20x_readSummary.RData")
  hide <- hi[grep("DE,TP",hi$evi)]
  hi <- hi[-grep("DE,TP",hi$evi)]
  hi$TE <- gsub(":\\S*","",hi$TE)
  hide$TE <- gsub(":\\S*","",hide$TE)
  
  load(paste0(DATADIR,"CallsForHeatMaps_gm12878_20xsc.RData"))
  table(res$call)
  res$id <- paste0("id",1:length(res))
  m.unspl.hi <- c()
  m.spl.hi <- c()
  
  win=4010
  bw=10
  for(TE in c("Alu","SVA","L1")){
    cat(TE,"\n")
    a1 <-res[res$TE==TE]
    ind <- resize(a1,win,"center")
    ind <- unlist(tile(ind,width = bw))
    
    gr2 <- hi[hi$TE%in%c(TE,"PolyA")] #,"PolyA"
    #if(TE=="L1"){
    #  gr1 <- gr1[gr1$mm<3]
    #  gr2 <- gr2[gr2$mm<3]
    #}
    
    ind$hi <- countOverlaps(ind,gr2,ignore.strand=T)
    ind$hide <- countOverlaps(ind,hide[hide$TE==TE],ignore.strand=T) # | hide$TE=="PolyA"
    
    x1 <- as.data.frame(matrix(ind$hide,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    y1 <- as.data.frame(matrix(ind$hi,nrow = length(a1),ncol = round(win/bw), byrow = T) )
    rownames(x1) <- rownames(y1) <- a1$id
    
    m.unspl.hi <- rbind(m.unspl.hi,y1)
    m.spl.hi <- rbind(m.spl.hi,x1)
    
    rm(x,y,x1,y1,ind,a1,gr1,gr2)
  }
  rm(wgs,hi,wgsde,hide)
  m.spl.hi <- m.spl.hi[match(res$id,rownames(m.spl.hi)),]
  m.unspl.hi <- m.unspl.hi[match(res$id,rownames(m.unspl.hi)),]
  
  save(res,m.spl.hi,m.unspl.hi,file=paste0(DATADIR,"CoveragePlotMatrices_20Xsc.RData"))
  
}


