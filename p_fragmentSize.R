

source("vars.R")
source("functions.R")


################################
if(F){
  myf<-function(path){
    x <-read.delim(gzfile(path),header=F)
    x <- data.table(x)
    x <- x[,total:=sum(V3),by=list(V1)]
    head(x)
    x <- as.data.frame(x)
    x <- unique(x[,c(1,4)])
    x 
    
  }
  
  x <- myf(paste0(RAWDIR,"gm12878_50xsc.ReadPairInsertSizeOriSummary.logs.gz"))
  y <- myf(paste0(RAWDIR,"gm12878_wgs_sc.ReadPairInsertSizeOriSummary.logs.gz"))
  
  pl <- data.frame(dis=seq(0,10,0.1))
  pl <- merge(pl,x,by.x="dis",by.y="V1",all.x=T)
  pl <- merge(pl,y,by.x="dis",by.y="V1",all.x=T)
  names(pl)[2:3] <- c("Hi-C","WGS")
  pl[is.na(pl)] <- 0
  names(pl)[2:3] <- c("Hi-C","WGS")
  #pl <- pl[pl$`Hi-C`+pl$WGS>0,]
  pl$`Hi-C` <- round(pl$`Hi-C`*100/sum(pl$`Hi-C`),2)
  pl$WGS <- round(pl$WGS*100/sum(pl$WGS),5)
  pl <- pl[1:90,]  
  pl <- melt(pl,measure.vars = names(pl)[2:3])
  pl$dis <-10^pl$dis
  rect1 <- data.frame (xmin=10, xmax=1000, ymin=0, ymax=Inf)
  rect2 <- data.frame (xmin=1000, xmax=Inf, ymin=0, ymax=Inf)
  pl$variable <- factor(pl$variable,levels=c("WGS","Hi-C"))
  pl <- pl[pl$variable]
  pl1 <- pl[(pl$value>0 & pl$variable=="Hi-C") | pl$variable=="WGS",]
  
  pdf(file=paste0(FIGDIR,"PF_InsertSizeDistribution.pdf"),width = 8,height = 4)
  ggplot(pl1, aes(x=dis, y=value))+geom_line(size=1.2,col="black")+
    facet_wrap(~variable,ncol=2,scales = "free_y")+
    ylab("% RPs")+xlab("effective insert sizes")+
    scale_color_brewer(type = "div",palette = "Dark2")+
    scale_y_continuous(expand = c(0,0))+theme_minimal()+
    geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray", alpha=0.3, inherit.aes = FALSE)+
    geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray40", alpha=0.3, inherit.aes = FALSE)+
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),expand=c(0,0))+
    theme(axis.title = element_text(size = 20,colour = "black"),
          axis.text = element_text(size = 20,colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 22,color="black",face = 3))
  
  dev.off()  
  
}


