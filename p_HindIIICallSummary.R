rm(list=ls())
source("vars.R")
source("functions.R")

pllist <- list()

###################### call summary
if(T){
  load(file=paste0(DATADIR,"CallsForHeatMaps_GM12878_HindIII_r1c.RData"))
  gm50 <- mycalls(paste0(RAWDIR,"GM12878_HindIII_r1c.candidate.insertions.bed"),filt = F)
  gm50 <- gm50[gm50$TE%in%c("Alu","L1","SVA")]
  
  e <- res[res$call%in%c("PacB Germline","PacB only")]
  e1 <- subsetByOverlaps(gm50,resize(e,101,"center"),ignore.strand=T)
  v <- as.character(e1$remark)
  table(v)
  v[grep("NoEnrichment",v)] <- "no RAM enrichment"
  v[grep("<2.5%RAMCov|<5%ClipCov|RAM=0|InsuffClipCov",v)] <- "coverage thresholds"
  v[grep("Unmap",v)] <- "bad reciprocal cluster"
  v[grep("Amb|<2TE-mapping|poly|indel|<5-refMapq|ClipAt|<30Map|Clonal",v)] <- "miscellaneous filters"
  v[grep("noTail",v)] <- "tail absent"
  table(v)
  v <- c(v,rep("<2 TE-mapping clip reads/RE proximity", (length(e)-length(v))))
  
  
  ## call summary
  pl2 <- as.data.frame(table(res[!res$call%in%c("PacB Germline","PacB only")]$call))
  names(pl2)[1:2] <- c("status","numbers")
  pl2$Label <- paste0(pl2$status,"\n(",pl2$numbers,")")
  p <- ggplot(pl2, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity",color="black")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq",palette = "Oranges",direction = 1)+
    theme_void() + 
    geom_text_repel(aes(x=2,y=cumsum(pl2$numbers) - pl2$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  pllist[["a"]] <- print(p)
  
    ## reason of missed calls
  pl1 <- as.data.frame(table(v))
  names(pl1)[1:2] <- c("status","numbers")
  pl1$Label <- paste0(pl1$status,"\n(",pl1$numbers,")")
  p <- ggplot(pl1, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity",color="black")+
    coord_polar("y", start=0)+scale_fill_brewer(type = "seq",palette = "Oranges",direction = 1)+
    theme_void() + 
    geom_text_repel(aes(x=2,y=cumsum(pl1$numbers) - pl1$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
    theme(legend.position = "none",axis.ticks.x = element_line(size = 2))
  pllist[["b"]] <- print(p)
  
  rm(e,e1,res,p,gm50,v)
}

### coverage
if(T){
  load(paste0(DATADIR,"Coverage.RData"))
  
  d <- zz[zz$profile%in%c("DpnII (1.4B)", "WGS (1.4B)","HindIII (1.8B)"),]
  d$profile <- factor(d$profile,levels=c("DpnII (1.4B)","HindIII (1.8B)", "WGS (1.4B)"))
  d <- d[order(d$profile,d$coverage),]
  d$profile <- as.character(d$profile)
  d$coverage <- as.numeric(d$coverage)
  d$totper <- as.numeric(d$totper)
  e <- reshape2::dcast(d,coverage~profile,value.var = "totper")
  e[,2:ncol(e)] <- apply(e[,2:ncol(e)],2,cumsum)
  e <- melt(e,measure.vars = names(e)[2:ncol(e)])
  names(e)[2] <- "profile"
  e$profile <- factor(e$profile,levels=c("DpnII (1.4B)","HindIII (1.8B)", "WGS (1.4B)"))
  rect1 <- data.frame (xmin=30, xmax=31, ymin=0, ymax=Inf)
  e$value <- 100 - e$value
  
  cols <- c("#3182BD","#62D469", "red"  )
  names(cols) <- c("DpnII (1.4B)","HindIII (1.8B)", "WGS (1.4B)")
  
  p <- ggplot(e, aes(x=(coverage), y=(value),col=profile)) + geom_line(size=1) + geom_point(size=1.1)+
    xlab("coverage (per 50bp bin)") + ylab("% of genome (cumulative)")+scale_color_manual(values = cols)+
    ylim(0,100)+
    scale_x_continuous(expand = c(0,0),breaks = c(0,25,50,75,100),labels = c(0,25,50,75,">100"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,100))+
    geom_vline(xintercept = 30,color="gray",linetype="dashed",size=1.1)+
    theme_minimal()+theme(axis.line = element_line(color="gray",size=0.8),
                          axis.ticks = element_line(colour = "gray"),
                          axis.title = element_text(size=16,color="black"),
                          axis.text = element_text(size=15,colour = "black"),
                          legend.text = element_text(size=14,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "bottom")
  pllist[["c"]] <-print(p)
  rm(d,p,rect1,zz)
}

### SPSN
if(T){
  load(paste0(DATADIR,"SPSN_dataframe.RData"))
  col=c("MboI_precision"="#bcbddc","MboI_recall"="#efedf5",
        "precision"="#a1d99b","recall"="#e5f5e0",
        "WGS_precision"="#de2d26","WGS_recall"="#fee0d2")
  
  myres3 <- c()
  for(f in c("Alu","L1", "SVA")){
    res1 <- res[res$profile%in%c("1.4B (MboI)", "WGS_MELT","WGS_MELT_Alt","1.8B (HindIII)"),]
    
    res1 <- res1[res1$TE==f,]
    res1 <- res1[,c("profile","recall","precision")]
    res1 <- melt(res1,measure.vars = names(res1)[2:3])
    res1$profile <- factor(res1$profile,levels=c("1.8B (HindIII)","1.4B (MboI)",  "WGS_MELT", "WGS_MELT_Alt"))
    
    res1$variable <- gsub("fdr","FDR",res1$variable)
    res1$variable <- factor(res1$variable,levels=rev(c("recall","precision")))
    res1$TE <- f
    myres3 <- rbind(myres3,res1)
    rm(res1)
  }
  myres3$variable <- as.character(myres3$variable)
  myres3[grep("MELT",myres3$profile),]$variable <- paste0("WGS_",myres3[grep("MELT",myres3$profile),]$variable)
  myres3[grep("MboI",myres3$profile),]$variable <- paste0("MboI_",myres3[grep("MboI",myres3$profile),]$variable)
  
  
  p <- ggplot(myres3,aes(x=profile, y=value,fill=variable,group=variable))+
    geom_bar(stat="identity",position='dodge',color="black",size=0.5)+
    facet_wrap(~TE,ncol=1)+scale_fill_manual(values = col)+
    theme_minimal()+theme(axis.text = element_text(size=14,colour = "black"),
                          axis.line.y = element_line(color = "gray"),
                          axis.ticks.y = element_line(color = "gray"),
                          axis.title = element_blank(),
                          legend.text = element_text(size=14,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = c(0.8,0.65),
                          strip.text = element_text(size=16,face = 3,color="black",hjust = 0),
                          strip.background = element_blank())
  pllist[["d"]] <- print(p)
  rm(p,res,myres3)
}


### combined plot
if(F){
  lay <- rbind(c(3,1),
               c(4,2),
               c(4,NA))
  
  pdf(paste0(FIGDIR,"CallSummary_missedCalls_GM12878HindIII.pdf"),width = 6,height = 8)
  grid.arrange(grobs = pllist,layout_matrix=lay)
  dev.off()
  
  
}


######################
######################