rm(list=ls())
source("functions.R")
source("vars.R")


if(F){
  d <- rbind(read.delim("C:/d/coverage/10xsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/20xsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/50xsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/wgssummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/100xsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/3xsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/HindIIIR1summarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/K562_WGSsummarizedCoverage.txt",header=F),
             read.delim("C:/d/coverage/K562summarizedCoverage.txt",header=F))
  names(d) <- c("profile","coverage","percent")
  d1 <- d[d$profile=="wgs",]
  plot(log2(d1$coverage+1),d1$percent)
  e <- d
  e$totreads <- as.numeric(e$coverage*e$percent)
  e <- unique(e)
  e <- e[,c(1,4)]
  e <- data.table(e)
  e <- e[,tot:=sum(totreads),by=list(profile)]
  e <- unique(e[,c(1,3)])
  e$tot <- round(e$tot/1e6,2)
  
  d$percent <- d$percent*100/61765409
  d$profile <- gsub("HindIIIR1","HindIII (1.8B)",d$profile)
  d$profile <- gsub("K562_WGS","WGS (1.1B)",d$profile)
  d$profile <- gsub("K562","DpnII (1.1B)",d$profile)
  d$profile <- gsub("hind3","HindIII (1.2B)",d$profile)
  d$profile <- gsub("nco1","NcoI (600M)",d$profile)
  d$profile <- gsub("100x","DpnII (5B)",d$profile)
  d$profile <- gsub("3x","DpnII (100M)",d$profile)
  d$profile <- gsub("20x","DpnII (600M)",d$profile)
  d$profile <- gsub("50x","DpnII (1.4B)",d$profile)
  d$profile <- gsub("10x","DpnII (300M)",d$profile)
  d$profile <- gsub("wgs","WGS (1.4B)",d$profile)
  
  d$coverage <- ifelse(d$coverage>101,101,d$coverage)
  d <- data.table(d)
  d <- d[,totper:=sum(percent),by=list(profile,coverage)]
  d <- as.data.frame(d)
  d <- unique(d[,c(1,2,4)])
  zz <- d
  save(zz,file=paste0(DATADIR,"Coverage.RData"))
}