rm(list=ls())
source("vars.R")
source("functions.R")

pllist <- list()
### Figure2

 ## spsn
if(T){
  load(paste0(DATADIR,"SPSN_dataframe.RData"))
  col=c("precision"="#bcbddc","recall"="#efedf5",
        "WGS_precision"="#de2d26","WGS_recall"="#fee0d2")
  
  
  #### barplot
  myres1 <- c()
  myres2 <- c()
  for(f in c("Alu","L1", "SVA")){
    res1 <- res[res$profile%in%c("1.4B (MboI)", "300M (MboI)", "5B (MboI)", 
                                 "600M (MboI)", "WGS_MELT","WGS_MELT_Alt"),]
    res1$profile <- gsub("[(]MboI[)]","",res1$profile)
    res1$profile <- gsub("100M","100M\nHiTEA\nHi-C",res1$profile)
    res1$profile <- gsub("300M","300M\nHiTEA\nHi-C",res1$profile)
    res1$profile <- gsub("600M","600M\nHiTEA\nHi-C",res1$profile)
    res1$profile <- gsub("1.4B","1.4B\nHiTEA\nHi-C",res1$profile)
    res1$profile <- gsub("5B","5B\nHiTEA\nHi-C",res1$profile)
    res1$profile <- gsub("WGS_MELT_Alt","1.4B\nMELT-GT\nWGS",res1$profile)
    res1$profile <- gsub("WGS_MELT","1.4B\nMELT-PASS\nWGS",res1$profile)
    
    res2 <- res[res$profile%in%c("1.4B (MboI)", "WGS_MELT","WGS_MELT_Alt"),]
    res2$profile <- gsub("1.4B [(]MboI[)]","1.4B\nHiTEA\nHi-C",res2$profile)
    res2$profile <- gsub("WGS_MELT_Alt","1.4B\nMELT-GT\nWGS",res2$profile)
    res2$profile <- gsub("WGS_MELT","1.4B\nMELT-PASS\nWGS",res2$profile)
    
    res1 <- res1[res1$TE==f,]
    res2 <- res2[res2$TE==f,]
    res1 <- res1[,c("profile","recall","precision")]
    res2 <- res2[,c("profile","recall","precision")]
    res1 <- melt(res1,measure.vars = names(res1)[2:3])
    res2 <- melt(res2,measure.vars = names(res2)[2:3])
    res1$profile <- factor(res1$profile,levels=c("300M\nHiTEA\nHi-C ", "600M\nHiTEA\nHi-C ", "1.4B\nHiTEA\nHi-C ", 
                                                 "5B\nHiTEA\nHi-C ", "1.4B\nMELT-PASS\nWGS", "1.4B\nMELT-GT\nWGS"))
    res2$profile <- factor(res2$profile,levels=c("1.4B\nHiTEA\nHi-C", "1.4B\nMELT-PASS\nWGS", "1.4B\nMELT-GT\nWGS"))
    
    res1$variable <- gsub("fdr","FDR",res1$variable)
    res1$variable <- factor(res1$variable,levels=rev(c("recall","precision")))
    res2$variable <- gsub("fdr","FDR",res2$variable)
    res2$variable <- factor(res2$variable,levels=rev(c("recall","precision")))
    
    #res2 <- res2[res2$profile!="1.2B\nMELT\nWGS",]
    res1 <- res1[-grep("MELT",res1$profile),]
    res1$TE <- f
    res2$TE <- f
    myres1 <- rbind(myres1,res1)
    myres2 <- rbind(myres2,res2)
    rm(res1,res2)
  }
  myres2$variable <- as.character(myres2$variable)
  myres2[grep("MELT",myres2$profile),]$variable <- paste0("WGS_",myres2[grep("MELT",myres2$profile),]$variable)
  #myres2$variable <- factor(myres2$variable,levels=c("precision","recall"))
  
  p <- ggplot(myres1,aes(x=profile, y=value,fill=variable,group=variable))+
    geom_bar(stat="identity",position='dodge',color="black",size=0.5)+
    facet_wrap(~TE,ncol=1)+scale_fill_manual(values = col)+
    theme_minimal()+theme(axis.text = element_text(size=12,colour = "black"),
                          axis.line.y = element_line(color = "gray"),
                          axis.ticks.y = element_line(color = "gray"),
                          axis.title = element_blank(),
                          legend.text = element_text(size=12,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "none",
                          strip.text = element_text(size=16,face = 3,color="white",hjust = 0),
                          strip.background = element_blank())
  pllist[["a"]] <- print(p)
  p <- ggplot(myres2,aes(x=profile, y=value,fill=variable,group=variable))+
    geom_bar(stat="identity",position='dodge',color="black",size=0.5)+
    facet_wrap(~TE,ncol=1)+scale_fill_manual(values = col)+
    theme_minimal()+theme(axis.text = element_text(size=12,colour = "black"),
                          axis.line.y = element_line(color = "gray"),
                          axis.ticks.y = element_line(color = "gray"),
                          axis.title = element_blank(),
                          legend.text = element_text(size=12,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = c(0.8,0.65),
                          strip.text = element_text(size=16,face = 3,color="black",hjust = 0),
                          strip.background = element_blank())
  pllist[["b"]] <- print(p)
  
  rm(p,res,f)
}


## cov
if(T){
  load(paste0(DATADIR,"Coverage.RData"))
  ## cumsum, plot1
  d <- zz[zz$profile%in%c("DpnII (100M)","DpnII (300M)", "DpnII (600M)", "DpnII (1.4B)", "WGS (1.4B)", 
                          "DpnII (5B)"),]
  d <- d[order(d$profile,d$coverage),]
  d$profile <- as.character(d$profile)
  d$coverage <- as.numeric(d$coverage)
  d$totper <- as.numeric(d$totper)
  e <- reshape2::dcast(d,coverage~profile,value.var = "totper")
  e[,2:ncol(e)] <- apply(e[,2:ncol(e)],2,cumsum)
  e <- melt(e,measure.vars = names(e)[2:ncol(e)])
  names(e)[2] <- "profile"
  e$profile <- factor(e$profile,levels=c("DpnII (100M)","DpnII (300M)","DpnII (600M)",
                                         "DpnII (1.4B)","DpnII (5B)", "WGS (1.4B)"))
  rect1 <- data.frame (xmin=30, xmax=31, ymin=0, ymax=Inf)
  e$value <- 100 - e$value
  
  cols <- c(brewer.pal(5,"Blues"),"red")
  names(cols) <- c("DpnII (100M)","DpnII (300M)","DpnII (600M)",
                   "DpnII (1.4B)","DpnII (5B)", "WGS (1.4B)")
  
  p <-ggplot(e, aes(x=(coverage), y=(value),col=profile)) + geom_line(size=1) + geom_point(size=0.8)+
    xlab("coverage (per 50bp bin)") + ylab("% of genome (cumulative)")+scale_color_manual(values = cols)+
    ylim(0,100)+
    scale_x_continuous(expand = c(0,0),breaks = c(0,25,50,75,100),labels = c(0,25,50,75,">100"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,100))+
    geom_vline(xintercept = 30,color="gray",linetype="dashed",size=0.8)+
    theme_minimal()+theme(axis.line = element_line(color="gray",size=0.8),
                          axis.ticks = element_line(colour = "gray"),
                          axis.title = element_text(size=13,color="black"),
                          axis.text = element_text(size=12,colour = "black"),
                          legend.text = element_text(size=12,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "bottom")
  pllist[["c"]] <- print(p)
  #geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray", alpha=0.3, inherit.aes = FALSE)
  rm(d,zz,rect1,p)
}



## combine
if(F){
  lay <- rbind(c(rep(3,6)),
               c(2,2,2,1,1,1))
  
  pdf(paste0(FIGDIR,"PF_Fig2.pdf"),width = 5,height = 8)
  grid.arrange(grobs = pllist,layout_matrix=lay,heights=c(3,6))
  dev.off()
  
  
}

