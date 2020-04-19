rm(list=ls())
source("functions.R")
source("vars.R")

### Paper figures
## where are pacb germline calls?
if(F){
  load("A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_50xsbpacBBHiT.RData")
  path="C:/d/gm12878_50xsb.candidate.insertions.bed"
  gm50 <- mycalls(path,filt = F)
  gm50 <- gm50[gm50$EVI>0]
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  
  ## Repeat annotation
  load("A:/work/hg38_annotations/Hg38_repeats.RData")
  o = findOverlaps(gm50,hg38,ignore.strand=T)
  gm50$rep <- NA
  gm50[queryHits(o)]$rep <- as.character(hg38[subjectHits(o)]$repFamily)
  o = findOverlaps(res,hg38,ignore.strand=T)
  res$rep <- NA
  res[queryHits(o)]$rep <- as.character(hg38[subjectHits(o)]$repFamily)
  
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
  
  
  ## calls overlapping to SVs in GM12878
  load("C:/d/gm12878_pacB/HG001_GRCh38_PBSV_Indels_usableSV.RData")
  sv <- sv[sv$type%in%c("INS")]
  #sv <- sv[width(sv)>50]
  res$sv <- countOverlaps(resize(res,101,"center"),sv,ignore.strand=T)
  res$sv <- ifelse(res$sv>1,1,res$sv)
  o <- findOverlaps(resize(res,101,"center"),sv,ignore.strand=T)
  res$SVtype <- ""
  res[queryHits(o)]$SVtype <- sv[subjectHits(o)]$V6
  #res$SVSize <- 0
  #res[queryHits(o)]$SVsize <- width(sv[subjectHits(o)])
  
  
  ## check if the discordent reads are enriched on both sides
  load(file="C:/d/WGS_MELTRAM.RData")
  wgs <- wgs[-grep("S|H",wgs$cigar)]
  p <- GRanges()
  for(f in c("Alu","L1","SVA")){
    cat(f,"\n")
    g <- wgs[wgs$te==f]
    r <- res[res$TE==f]
    r$up <- countOverlaps(flank(r,200),g[strand(g)=="+"],ignore.strand=T) 
    r$down <- countOverlaps(flank(r,200,start = F),g[strand(g)=="-"],ignore.strand=T) 
    p <- c(p,r)
    rm(g,r)
  }
  res <- p
  rm(p)
  res$test <- ifelse(res$down/(res$up+res$down)<0.01 | res$up/(res$up+res$down) <0.01|(res$up+res$down)<5,1,0)
  sum(res$test)
  table(res$test,res$melt_pass)
  
  v <- c()
  n <- c()
  for(f in levels(factor(res$call))){
    v <- append(v,sum(table(as.character(res[res$call==f]$rep)))*100/length(res[res$call==f]))
    n <- append(n,f)
  }
  n <- c("HiTEA only\n(Overlap Germline calls)", "HiTEA\n(Overlap PcBio calls)","PacBio only\n(Overlap Germline calls)",
         "PacBio only","HiTEA only")
  pl <- data.frame(calls=n,percent=v)
  pl$calls <- factor(pl$calls,levels=rev(c("HiTEA only\n(Overlap Germline calls)","HiTEA only", "HiTEA\n(Overlap PcBio calls)","PacBio only","PacBio only\n(Overlap Germline calls)")))
  
  ### % calls in the repeat region
  #pdf('A:/work/RT in HiC/avg_int/CallsMadeInRepeatRegions_byCategory.pdf',width = 6,height = 5)
  ggplot(pl,aes(x=calls,y=percent,fill=calls))+geom_bar(stat="identity") +coord_flip()+
    scale_fill_brewer(palette="Dark2")+xlab("")+ylab("% of calls \nin the repeat regions")+
    geom_text(aes(label=paste0(round(percent,2),"%")), hjust=1,vjust=0.5,size=6) +
    theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                     axis.text = element_text(size=16,colour = "black"),
                     legend.text = element_text(size=14,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_blank())
  #dev.off()
  
  res$rep <- ifelse(is.na(res$rep),"",as.character(res$rep))
  res$sameR <- 0
  res[res$TE==res$rep]$sameR <- 1
  r <- res[res$sameR==1]
  table(r$call)*100/table(res$call)
  #r <- r[r$call!="unique"]
  length(r[r$melt_pass>0])*100/length(res[res$melt_pass>0])
  length(r[r$melt_gt>0])*100/length(res[res$melt_gt>0])
  length(r[r$call%in%c("HiTEA+PacB")])*100/length(res[res$call%in%c("HiTEA+PacB")])
  length(r[!r$call%in%c("Germline","unique")])
  
  v <- c(length(r[r$melt_pass>0]),length(r[r$melt_gt>0]),length(r[r$call%in%c("HiTEA+PacB")]))
  pl <- data.frame(value=v,
                   caller=c("MELT-PASS","MELT-GT","HiTEA"))
  pl$label <- paste0(pl$value,"/",length(r[!r$call%in%c("Germline","unique")]))
  
  #pdf('A:/work/RT in HiC/avg_int/PF_CallsWithinSameRef copy.pdf',width = 6,height = 4)
  ggplot(pl,aes(x=caller,y=value))+geom_bar(stat="identity",fill="lightblue") +
    ylab("number of true-positive insertions\n within same reference copy")+xlab("")+
    coord_flip()+
    geom_text(aes(label=label), hjust=0.2,vjust=0.5,size=6) +ylim(0,85)+
    theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                     axis.text = element_text(size=16,colour = "black"),
                     legend.text = element_text(size=14,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_blank())
  #dev.off()
  
  
  
  
  ## %calls within the SV containing locus
  if(F){
    pl2 <- as.data.frame(table(res[res$sv>0]$call)*100/table(res$call))
    pl2$Var1 <- factor(pl2$Var1,levels=rev(c("Germline", "unique", "HiTEA+PacB", "PacB only", "PacB Germline")))
    r <- res[res$sv>0]
    length(r[r$melt_pass>0])*100/length(res[res$melt_pass>0])
    length(r[r$melt_gt>0])*100/length(res[res$melt_gt>0])
    length(r[r$call%in%c("HiTEA+PacB")])*100/length(res[res$call%in%c("HiTEA+PacB")])
    length(r[!r$call%in%c("Germline","unique")])
    
    pdf('A:/work/RT in HiC/avg_int/CallsMadeInSVRegions_byCategory.pdf',width = 6,height = 5)
    ggplot(pl2,aes(x=Var1,y=Freq,fill=Var1))+geom_bar(stat="identity") +coord_flip()+
      scale_fill_brewer(palette="Dark2")+xlab("")+ylab("% of calls \nalong SV locus")+
      geom_text(aes(label=paste0(round(Freq,2),"%")), hjust=0,vjust=0.5,size=6) +ylim(0,45)+
      theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                       axis.text = element_text(size=16,colour = "black"),
                       legend.text = element_text(size=14,colour = "black"),
                       legend.title = element_blank(),
                       legend.position = "none",
                       strip.text = element_blank())
    dev.off()
    
    res$lowmelt <- countOverlaps(resize(res,101,"center"),m[m$Evi!="PASS"],ignore.strand=T)
    res$lowmelt <- ifelse(res$lowmelt>1,1,res$lowmelt)
    table(res$call,res$lowmelt)
    
    r <- res[res$TE==res$rep]
    table(r$call,r$melt_gt)
    table(r$call,r$test)
  }
  
  
  g <- res[res$call=="HiTEA+PacB" & res$melt_pass==0]
  e <- res[res$call%in%c("PacB Germline","PacB only") & res$melt_pass>0]
  g <- res[res$call=="HiTEA+PacB" & res$melt_gt==0]
  e <- res[res$call%in%c("PacB Germline","PacB only") & res$melt_gt>0]
  
  
  cf <- rbind(data.frame(TE=g$TE, category="missed by MELT (WGS)"),
              data.frame(TE=e$TE, category="missed by HiTEA (Hi-C)"))
  cf$TE <- factor(cf$TE,levels=rev(c("Alu","L1","SVA")))
  ## 
  pdf('A:/work/RT in HiC/avg_int/PF_MissRates.pdf',width = 8,height = 4)
  ggplot(cf,aes(x=TE,fill=TE))+geom_bar(stat="count")+
    facet_wrap(~category,ncol=2,scales = "free")+
    scale_fill_brewer(palette="Dark2")+xlab("")+
    geom_text(stat='count', aes(label=..count..), vjust=0.5,hjust=0.5,size=6) +
    theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                     axis.text = element_text(size=18,colour = "black"),
                     legend.text = element_text(size=14,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_text(size=18,colour = "black"),
                     strip.background = element_blank())
  dev.off()
  
  
  ## what are the MELT missed calls?
  v <- c(length(g[g$rep==g$TE & g$test==0]),length(g[g$test>0 & g$sv==0]),
         length(g[g$sv>0]))
  v <- c(v,length(g)-sum(v))
  
  pl1 <- data.frame(numbers=v,status=c("repeat olap","one sided RAM","proximal to SV","unknown"))
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  
  pdf('A:/work/RT in HiC/avg_int/PF_MELTMissed_call_Remarks.pdf',width = 4,height = 4)
  ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq",direction = -1)+
    theme_void() +
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  dev.off()  
  
  
  
  ## which filters the missed calls were removed by?
  gm50 <- mycalls(path,filt = F)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  e1 <- subsetByOverlaps(gm50,resize(e,101,"center"),ignore.strand=T)
  v <- as.character(e1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("5xRAMCov|10xClipCov|10xClipTECov",v)] <- "coverage thresholds"
  v[grep("poorTECluster;",v)] <- "bad reciprocal cluster"
  v[grep("PoorTEMapScore|fuzzy|doubleClip|badPair;",v)] <- "miscellaneous filters"
  v <- c(v,rep("<2 TE-mapping clip reads", (length(e)-length(v))))
  
  pl1 <- as.data.frame(table(v))
  names(pl1)[1:2] <- c("status","numbers")
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  
  pdf('A:/work/RT in HiC/avg_int/PF_HiTEAMissed_call_Remarks.pdf',width = 4,height = 4)
  ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq")+
    theme_void() +
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  dev.off()  
  
  pie(table(v))
  barplot(table(as.character(g$TE)),las=1,main="missed by MELT (WGS)")
  barplot(table(e$TE),las=1,main="missed by HiTEA (Hi-C)")
  
  
}

### Paper figures
######### HindIII_r1 call summary as supple fig
if(T){
  load(file="A:/work/RT in HiC/avg_int/CallsForHeatMaps_Gm12878_HindIII_r1pacBBHiT.RData")
  path="C:/d/GM12878_HindIII_r1.candidate.insertions.bed"
  gm50 <- mycalls(path,filt = F)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  pie(table(res[!res$call%in%c("PacB Germline","PacB only")]$call))
  
  e <- res[res$call%in%c("PacB Germline","PacB only")]
  e1 <- subsetByOverlaps(gm50,resize(e,101,"center"),ignore.strand=T)
  v <- as.character(e1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("5xRAMCov|10xClipCov|10xClipTECov",v)] <- "failed due to coverage thresholds"
  v[grep("poorTECluster;",v)] <- "bad reciprocal cluster"
  v[grep("PoorTEMapScore|fuzzy|doubleClip|badPair|bgOlap;",v)] <- "miscellaneous filters"
  v[grep("poorRecClust;",v)] <- "bad reciprocal clusters"
  v <- c(v,rep("<2 TE-mapping clip reads", (length(e)-length(v))))
  
  pl1 <- as.data.frame(table(v))
  names(pl1)[1:2] <- c("status","numbers")
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  
  pdf('A:/work/RT in HiC/avg_int/PF_HiTEAMissed_call_Remarks_HindIII_r1.pdf',width = 4,height = 4)
  ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq")+
    theme_void() +
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  dev.off()  
  
  
  pl1 <- as.data.frame(table(res[!res$call%in%c("PacB Germline","PacB only")]$call))
  names(pl1)[1:2] <- c("status","numbers")
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  
  pdf('A:/work/RT in HiC/avg_int/PF_HiTEAMissed_HindIII_r1_calls.pdf',width = 4,height = 4)
  ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq")+
    theme_void() +
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  dev.off()  
  
  
  pie(table(v))
  
}