rm(list=ls())
source("vars.R")
source("functions.R")


### Coverage
if(F){
  load(paste0(DATADIR,"Coverage_CallsMissedbyMELTHiTEA_GM1287850x.RData"))
  
  pdf(paste0(FIGDIR,"PF_MELT_HiTEAMissdCallCov_boxplot.pdf"),width = 5,height = 4)
  ggplot(z,aes(x=as.factor(def),y=change))+geom_violin()+ylim(-4,4)+
    xlab("")+ylab("coverage ratio (Hi-C/WGS), log2")+
    stat_compare_means(size=5)+
    geom_hline(yintercept = 0)+
    geom_boxplot(width=0.5,fill="#deebf7")+theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                                                              axis.text = element_text(size=16,colour = "black"),
                                                              legend.text = element_text(size=14,colour = "black"),
                                                              legend.title = element_blank(),
                                                              legend.position = "none",
                                                              strip.text = element_blank())
  dev.off()
  
}


### reasons why calls are missed
if(F){
  load(paste0(DATADIR,"CallsMissedAnnotations_GM1287850x.RData"))
  
  ### % calls in the repeat region
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
  
  ggplot(pl,aes(x=calls,y=percent,fill=calls))+geom_bar(stat="identity") +coord_flip()+
    scale_fill_brewer(palette="Dark2")+xlab("")+ylab("% of calls \nin the repeat regions")+
    geom_text(aes(label=paste0(round(percent,2),"%")), hjust=1,vjust=0.5,size=6) +
    theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                     axis.text = element_text(size=16,colour = "black"),
                     legend.text = element_text(size=14,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_blank())
  
  
  
  #### TE insertions within the same reference copy
  res$rep <- ifelse(is.na(res$rep),"",as.character(res$rep))
  res$sameR <- 0
  res[res$TE==res$rep]$sameR <- 1
  r <- res[res$sameR==1]
  table(r$call)*100/table(res$call)
  length(r[r$melt_pass>0])*100/length(res[res$melt_pass>0])
  length(r[r$melt_gt>0])*100/length(res[res$melt_gt>0])
  length(r[r$call%in%c("HiTEA+PacB")])*100/length(res[res$call%in%c("HiTEA+PacB")])
  length(r[!r$call%in%c("Germline","unique")])
  
  v <- c(length(r[r$melt_pass>0]),length(r[r$melt_gt>0]),length(r[r$call%in%c("HiTEA+PacB")]))
  pl <- data.frame(value=v,caller=c("MELT-PASS","MELT-GT","HiTEA"))
  pl$label <- paste0(pl$value,"/",length(r[!r$call%in%c("Germline","unique")]))
  
  cols =c("MELT-PASS"="#de2d26","MELT-GT"="#de2d26","HiTEA"="#bcbddc")
  
  pdf(paste0(FIGDIR,'PF_CallsWithinSameRef copy.pdf'),width = 6,height = 4)
  ggplot(pl,aes(x=caller,y=value,fill=caller))+geom_bar(stat="identity") +
    ylab("number of true-positive insertions\n within same reference copy")+xlab("")+
    coord_flip()+scale_fill_manual(values = cols)+
    geom_text(aes(label=label), hjust=0.2,vjust=0.5,size=6) +ylim(0,85)+
    theme_bw()+theme(axis.title = element_text(size=18,color="black"),
                     axis.text = element_text(size=16,colour = "black"),
                     legend.text = element_text(size=14,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_blank())
  dev.off()
  
  
  
  ### Numbers of calls missed by MELT or HiTEA 
  r <- res[res$TE==res$rep]
  table(r$call,r$melt_gt)
  table(r$call,r$RAMtest)
  
  #g <- res[res$call=="HiTEA+PacB" & res$melt_pass==0]
  #e <- res[res$call%in%c("PacB Germline","PacB only") & res$melt_pass>0]
  g <- res[res$call=="HiTEA+PacB" & res$melt_gt==0]
  e <- res[res$call%in%c("PacB Germline","PacB only") & res$melt_gt>0]

  cf <- rbind(data.frame(TE=g$TE, category="missed by MELT (WGS)"),
              data.frame(TE=e$TE, category="missed by HiTEA (Hi-C)"))
  cf$TE <- factor(cf$TE,levels=rev(c("Alu","L1","SVA")))
  cols <- c("Alu"="#2c7fb8","SVA"="#edf8b1","L1"="#7fcdbb")
  
  pllist <- list()
  
  ## 
  p <- ggplot(cf,aes(x=TE,fill=TE))+geom_bar(stat="count",color="black")+
    facet_wrap(~category,ncol=2,scales = "free")+
    scale_fill_manual(values = cols)+xlab("")+
    geom_text(stat='count', aes(label=..count..), vjust=0.5,hjust=0.5,size=6) +
    theme_bw()+theme(axis.title = element_text(size=16,color="black"),
                     axis.text = element_text(size=15,colour = "black"),
                     legend.text = element_text(size=15,colour = "black"),
                     legend.title = element_blank(),
                     legend.position = "none",
                     strip.text = element_text(size=16,colour = "black"),
                     strip.background = element_blank())
  pllist[["counts"]] <- print(p)
  pdf(paste0(FIGDIR,'GM12878_MissedCallsBarplot.pdf'),width = 8,height = 4)
  pllist[["counts"]]
  dev.off()
  
  
  ## what are the MELT missed calls?
  v <- c(length(g[g$rep==g$TE & g$RAMtest==0]),length(g[g$RAMtest>0]))
  v <- c(v,length(g)-sum(v))
  pl1 <- data.frame(numbers=v,status=c("repeat olap","one sided RAM","unknown"))
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  pl1$plot <- 1
  
  p <- ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity",color="black")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq",palette = "Oranges",direction = 1)+
    theme_void() + 
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  pllist[["a"]] <- print(p)
  
  gm50 <- mycalls(paste0(RAWDIR,"gm12878_50xsc.candidate.insertions.bed"),filt = F)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  e1 <- subsetByOverlaps(gm50,resize(e,101,"center"),ignore.strand=T)
  v <- as.character(e1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("<2.5%RAMCov|<5%ClipCov|RAM=0|InsuffClipCov",v)] <- "coverage thresholds"
  v[grep("Unmap",v)] <- "bad reciprocal cluster"
  v[grep("noTail",v)] <- "tail absent"
  v[grep("Amb|<2TE-mapping|poly|indel|<5-refMapq",v)] <- "miscellaneous filters"
  table(v)
  v <- c(v,rep("<2 TE-mapping clip reads/RE proximity", (length(e)-length(v))))
  pl2 <- as.data.frame(table(v))
  names(pl2)[1:2] <- c("status","numbers")
  pl2$Label <- paste0(pl2$status,"\n(",pl2$numbers,")")
  pl2$plot <- 2
  pl2 <- pl2[,c("numbers","status","Label","plot")]

  p <- ggplot(pl2, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity",color="black")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq",palette = "Oranges",direction = 1)+
    theme_void() + 
    geom_text_repel(aes(x=2,y=cumsum(pl2$numbers) - pl2$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  pllist[["b"]] <- print(p)
  
  lay <- rbind(c(1,1),
               c(2,3))
  
  pdf(paste0(FIGDIR,"CallSummary_missedCalls_GM12878MboI.pdf"),width = 8,height = 6)
  grid.arrange(grobs = pllist,layout_matrix=lay)
  dev.off()
  
  
  
}



#3 missed call: numbers
if(F){
  load(paste0(DATADIR,"Coverage_CallsMissedbyMELTHiTEA_GM1287850x.RData"))
  
  gm50 <- mycalls(paste0(RAWDIR,"gm12878_50xsc.candidate.insertions.bed"),filt = F)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  
  table(res$call)
  e <- res[res$call=="PacB Germline"]
  f <- res[res$call=="PacB only"]
  
  e1 <- subsetByOverlaps(gm50,resize(e,101,"center"),ignore.strand=T)
  v <- as.character(e1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("<2.5%RAMCov|<5%ClipCov|RAM=0|InsuffClipCov",v)] <- "coverage thresholds"
  v[grep("Unmap",v)] <- "bad reciprocal cluster"
  v[grep("noTail",v)] <- "tail absent"
  v[grep("Amb|<2TE-mapping|poly|indel|<5-refMapq",v)] <- "miscellaneous filters"
  table(v)
  v <- c(v,rep("<2 TE-mapping clip reads/RE proximity", (length(e)-length(v))))
  table(v)
  
  
  f1 <- subsetByOverlaps(gm50,resize(f,101,"center"),ignore.strand=T)
  v <- as.character(f1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("<2.5%RAMCov|<5%ClipCov|RAM=0|InsuffClipCov",v)] <- "coverage thresholds"
  v[grep("Unmap",v)] <- "bad reciprocal cluster"
  v[grep("noTail",v)] <- "tail absent"
  v[grep("Amb|<2TE-mapping|poly|indel|<5-refMapq",v)] <- "miscellaneous filters"
  table(v)
  v <- c(v,rep("<2 TE-mapping clip reads/RE proximity", (length(f)-length(v))))
  table(v)
  
  
  
  load(paste0(DATADIR,"CallsMissedAnnotations_GM1287850x.RData"))
  r <- res[res$call=="unique"]
  sum(r$melt_gt>0)
  sum(r$melt_pass>0)
  r <- mycalls.germlineAnnotations_gno(r)
  
}
