rm(list=ls())
source("functions.R")
source("vars.R")

if(F){
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
  
  
   s <-ggplot(e, aes(x=(coverage), y=(value),col=profile)) + geom_line(size=1) + geom_point(size=1.75)+
    xlab("coverage (per 50bp bin)") + ylab("% of genome (cumulative)")+scale_color_manual(values = cols)+
    ylim(0,100)+
    scale_x_continuous(expand = c(0,0),breaks = c(0,25,50,75,100),labels = c(0,25,50,75,">100"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,100))+
    geom_vline(xintercept = 30,color="gray",linetype="dashed",size=1.4)+
     theme_minimal()+theme(axis.line = element_line(color="gray",size=0.8),
                          axis.ticks = element_line(colour = "gray"),
                          axis.title = element_text(size=16,color="black"),
                          axis.text = element_text(size=15,colour = "black"),
                          legend.text = element_text(size=14,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "bottom")
    #geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray", alpha=0.3, inherit.aes = FALSE)
  
  pdf(paste0(FIGDIR,"50bpBin_coverage_cumsum.pdf"),width = 6,height = 5)
  s 
  dev.off()
  
  
  
  ## cumsum, plot1
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
  
  s <- ggplot(e, aes(x=(coverage), y=(value),col=profile)) + geom_line(size=1) + geom_point(size=1.75)+
    xlab("coverage (per 50bp bin)") + ylab("% of genome (cumulative)")+scale_color_manual(values = cols)+
    ylim(0,100)+
    scale_x_continuous(expand = c(0,0),breaks = c(0,25,50,75,100),labels = c(0,25,50,75,">100"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,100))+
    geom_vline(xintercept = 30,color="gray",linetype="dashed",size=1.4)+
    theme_minimal()+theme(axis.line = element_line(color="gray",size=0.8),
                          axis.ticks = element_line(colour = "gray"),
                          axis.title = element_text(size=16,color="black"),
                          axis.text = element_text(size=15,colour = "black"),
                          legend.text = element_text(size=14,colour = "black"),
                          legend.title = element_blank(),
                          legend.position = "bottom")
  pdf(paste0(FIGDIR,"50bpBin_coverage_cumsum_Hind3.pdf"),width = 6,height = 5)
  s 
  dev.off()
  
}

