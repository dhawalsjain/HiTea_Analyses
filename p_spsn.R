rm(list=ls())
source("functions.R")
source("vars.R")

if(F){
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
  
  p1 <- ggplot(myres1,aes(x=profile, y=value,fill=variable,group=variable))+
    geom_bar(stat="identity",position='dodge',color="black",size=0.5)+
    facet_wrap(~TE,ncol=1)+scale_fill_manual(values = col)+
    theme_minimal()+theme(axis.text = element_text(size=14,colour = "black"),
                          axis.line.y = element_line(color = "gray"),
                          axis.ticks.y = element_line(color = "gray"),
                          axis.title = element_blank(),
                          legend.text = element_text(size=14,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "none",
                          strip.text = element_text(size=16,face = 3,color="white",hjust = 0),
                          strip.background = element_blank())
  p2 <- ggplot(myres2,aes(x=profile, y=value,fill=variable,group=variable))+
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
  
  
  
  
  pdf(paste0(FIGDIR,"SPSN_barplot_v3.pdf"),width = 6,height = 7)
  grid.arrange(p2,p1,nrow=1, widths=c(3,4))
  dev.off()
  
  ## HindIII
  col=c("precision"="#a1d99b","recall"="#e5f5e0",
        "MboI_precision"="#bcbddc","MboI_recall"="#efedf5",
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
  
  
  p3 <- ggplot(myres3,aes(x=profile, y=value,fill=variable,group=variable))+
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
  
  pdf(paste0(FIGDIR,"SPSN_barplot_HindIII_v2.pdf"),width = 5,height = 6)
  p3
  dev.off()
  
  
  ## 600M reads
  myres4 <- c()
  for(f in c("Alu","L1", "SVA")){
    res4 <- res[res$profile%in%c("600M (MboI)", "WGS_MELT_600M","WGS_MELT_600M_Alt"),]
    res4$profile <- gsub("600M [(]MboI[)]","600M\nHiTEA\nHi-C",res4$profile)
    res4$profile <- gsub("WGS_MELT_600M_Alt","600M\nMELT-GT\nWGS",res4$profile)
    res4$profile <- gsub("WGS_MELT_600M","600M\nMELT-PASS\nWGS",res4$profile)
    
    res4 <- res4[res4$TE==f,]
    res4 <- res4[,c("profile","recall","precision")]
    res4 <- melt(res4,measure.vars = names(res4)[2:3])
    res4$profile <- factor(res4$profile,levels=c("600M\nHiTEA\nHi-C", "600M\nMELT-PASS\nWGS", "600M\nMELT-GT\nWGS"))
    
    res4$variable <- gsub("fdr","FDR",res4$variable)
    res4$variable <- factor(res4$variable,levels=rev(c("recall","precision")))
    
    res4$TE <- f
    myres4 <- rbind(myres4,res4)
    rm(res4)
  }
  myres4$variable <- as.character(myres4$variable)
  myres4[grep("MELT",myres4$profile),]$variable <- paste0("WGS_",myres4[grep("MELT",myres4$profile),]$variable)
  
  col=c("precision"="#bcbddc","recall"="#efedf5",
        "WGS_precision"="#de2d26","WGS_recall"="#fee0d2")
  
  
  p4 <- ggplot(myres4,aes(x=profile, y=value,fill=variable,group=variable))+
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
  
  pdf(paste0(FIGDIR,"SPSN_barplot_GM12878_600M.pdf"),width = 5,height = 6)
  p4
  dev.off()
  
  
}



if(F){
  d <- read.delim(paste0(DATADIR,"pacb_references_hg38_RE.txt"),header=F,sep = "\t",stringsAsFactors = F)
  d$dist <- abs(d$V2-d$V7)
  
  ggplot(d,aes(x=log10(dist+1),col=V6))+stat_ecdf(geom="point",size=2)+
    xlim(0,5)
  
  
}