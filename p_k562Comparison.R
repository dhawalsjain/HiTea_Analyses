
rm(list=ls())
source("functions.R")
source("vars.R")

#####################################################################################3

if(F){
  
  require(UpSetR)
  cf <- c()
  
  #### K562
  ht <- mycalls(path = paste0(RAWDIR,"K562_sc.candidate.insertions.bed"),filt = T)
  #ht <- mycalls(path = "C:/d/K562_sb.candidate.insertions.bed",filt = T)
  
  n<- read.delim("A:/work/RT in HiC/MELT/k562_5m/Insertions.txt",comment.char = "#",header=F)
  n$V5 <- gsub("<INS:ME:","",n$V5)
  n$V5 <- gsub(">","",n$V5)
  n$V5 = gsub("ALU","Alu",n$V5)
  n$V5 = gsub("LINE1","L1",n$V5)
  n =  with(n,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  n = n[n$Evi=="PASS"]
  
  ngt<- read.delim("A:/work/RT in HiC/MELT/k562_5m/Insertions.txt",comment.char = "#",header=F)
  ngt[grep("0/1|1/1|1/0",ngt$V10),]$V7 <- "PASS"
  ngt$V5 <- gsub("<INS:ME:","",ngt$V5)
  ngt$V5 <- gsub(">","",ngt$V5)
  ngt$V5 = gsub("ALU","Alu",ngt$V5)
  ngt$V5 = gsub("LINE1","L1",ngt$V5)
  ngt =  with(ngt,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  ngt = ngt[ngt$Evi=="PASS"]
  
  for(f in c("Alu","L1","SVA")){
    Alu <- c(ht[ht$TE==f,0],n[n$TE==f,0],ngt[ngt$TE==f,0])
    Alu <- resize(reduce(resize(Alu,101,"center"),ignore.strand=T),1,"center")
    
    Alu$ht <- countOverlaps(resize(Alu,101,"center"),ht[ht$TE==f],ignore.strand=T)
    Alu$melt <- countOverlaps(resize(Alu,101,"center"),n[n$TE==f],ignore.strand=T)
    Alu$meltgt <- countOverlaps(resize(Alu,101,"center"),ngt[ngt$TE==f],ignore.strand=T)
    Alu <- mycalls.germlineAnnotations(Alu)
    #Alu$GLp <- ifelse(Alu$GLp=="Somatic",0,1)
    #Alu$GLp <- ifelse( Alu$melt>0 | Alu$meltgt>0,0,Alu$GLp)
    #Alu$GLp <- ifelse(Alu$GLp==1,as.character("Germline"),as.character("Somatic"))
    cf <- rbind(cf, data.frame(TE=f,HiTEA=Alu$ht,PASS=Alu$melt,GT=Alu$meltgt, Pop=Alu$GLp,cellLine="K562"))
    rm(Alu)
  }
  rm(n,ngt,ht)  
  
  ############ GM12878
  ht <- mycalls(path = paste0(RAWDIR,"gm12878_50xsc.candidate.insertions.bed"),filt = T)
  
  n<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  n$V5 <- gsub("<INS:ME:","",n$V5)
  n$V5 <- gsub(">","",n$V5)
  n$V5 = gsub("ALU","Alu",n$V5)
  n$V5 = gsub("LINE1","L1",n$V5)
  n =  with(n,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  n = n[n$Evi=="PASS"]
  
  ngt<- read.delim("A:/work/RT in HiC/melt/gm12878_5m/Insertions.txt",comment.char = "#",header=F)
  ngt[grep("0/1|1/1|1/0",ngt$V10),]$V7 <- "PASS"
  ngt$V5 <- gsub("<INS:ME:","",ngt$V5)
  ngt$V5 <- gsub(">","",ngt$V5)
  ngt$V5 = gsub("ALU","Alu",ngt$V5)
  ngt$V5 = gsub("LINE1","L1",ngt$V5)
  ngt =  with(ngt,GRanges(V1,IRanges(V2,V2),"+",TE=V5,Evi=V7,Descr=V8,GL="germline"))
  ngt = ngt[ngt$Evi=="PASS"]
  
  for(f in c("Alu","L1","SVA")){
    Alu <- c(ht[ht$TE==f,0],n[n$TE==f,0],ngt[ngt$TE==f,0])
    Alu <- resize(reduce(resize(Alu,101,"center"),ignore.strand=T),1,"center")
    
    Alu$ht <- countOverlaps(resize(Alu,101,"center"),ht[ht$TE==f],ignore.strand=T)
    Alu$melt <- countOverlaps(resize(Alu,101,"center"),n[n$TE==f],ignore.strand=T)
    Alu$meltgt <- countOverlaps(resize(Alu,101,"center"),ngt[ngt$TE==f],ignore.strand=T)
    Alu <- mycalls.germlineAnnotations(Alu)
    #Alu$GLp <- ifelse(Alu$GLp=="Somatic",0,1)
    #Alu$GLp <- ifelse( Alu$melt>0 | Alu$meltgt>0,0,Alu$GLp)
    #Alu$GLp <- ifelse(Alu$GLp==1,as.character("Germline"),as.character("Somatic"))
    cf <- rbind(cf, data.frame(TE=f,HiTEA=Alu$ht,PASS=Alu$melt,GT=Alu$meltgt, Pop=Alu$GLp, cellLine="GM12878"))
    rm(Alu)
  }
  rm(ht,n,ngt)
  
  
  ############## plot
  
  names(cf)[3:4]<- c("MELT-PASS","MELT-GT")
  
  
  
  plot_l <- list()
  for(cl in c("GM12878","K562")){
    for(f in c("Alu","L1","SVA")){
      m <- cf[cf$TE==f & cf$cellLine==cl,]
      i<-0
      mylist<-list()
      vectorUniqueValue <- c("Somatic","Germline")
      while ( length(vectorUniqueValue)>0 ){
        i<-i+1
        mylist[[ i ]]<-list(query = elements, params = list("Pop",as.character(vectorUniqueValue)), 
                            active = T,query.name=as.character(vectorUniqueValue)[1])
        vectorUniqueValue<-vectorUniqueValue[-1]
      }
      upset(m, sets = names(m)[2:4], mb.ratio = c(0.5, 0.5), order.by = c("freq"),
            text.scale=rep(1.5,6),point.size=5,matrix.color = "#68B343",
            queries = mylist,query.legend="bottom",color.pal = 1)
      rm(i,mylist,vectorUniqueValue)
      grid.edit('arrange', name=paste0(f,"_",cl))
      vp <- grid.grab()
      plot_l[[paste0(f,"_",cl)]] <- vp
    }
  }
  pdf(paste0(FIGDIR,"GM12878_K562_MELTCallComparisons.pdf"),width = 15,height = 7)
  grid.arrange(grobs=plot_l,nrow=2)
  dev.off()
  
  
  sum(cf$cellLine=="GM12878" & cf$HiTEA>0 )
  sum(cf$cellLine=="GM12878" & cf$HiTEA>0 & (cf$`MELT-GT`>0 | cf$Pop=="Germline"))
  
}




