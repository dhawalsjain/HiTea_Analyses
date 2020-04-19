rm(list=ls())
source("functions.R")
source("vars.R")
makematrix <- function(a){
  a = a[,100:300]
  a <- as.matrix(a)
  a[a>8] = 8
  #a <-log2(a+1)
  dim(a) = dim(a)
  attr(a, "upstream_index") = 1:101
  attr(a, "target_index") = 0
  attr(a, "downstream_index") = 101:200
  attr(a, "extend") = c(1000,1000)  # it must be a vector of length two
  class(a) = c("normalizedMatrix", "matrix")
  attr(a, "signal_name") = "Spt6"
  attr(a, "target_name") = "Insert"
  attr(a,"target_is_single_point") = TRUE
  return(a)
}

########### 50X #############################
if(F){
  library(EnrichedHeatmap)
  library(circlize)
  
  load(paste0(DATADIR,"CoveragePlotMatrices_50Xsc.RData"))
  table(res$call)
  m.spl.hi <- makematrix(m.spl.hi)
  m.spl.wgs <- makematrix(m.spl.wgs)
  m.unspl.hi <- makematrix(m.unspl.hi)
  m.unspl.wgs <- makematrix(m.unspl.wgs)
  maxval <- c(max(m.spl.hi),max(m.spl.wgs),max(m.unspl.hi),max(m.unspl.wgs))
  col = colorRamp2(c(0, max(maxval)), c("white","black"))
  
  res <- as.data.frame(res)
  res$call <- gsub("Germline","HiTEA (Germline)",res$call)
  res$call <- gsub("PacB HiTEA [(]Germline[)]","PacB (Germline)",res$call)
  table(res$call)
  
  partition = res$call
  partition <- factor(partition,levels=c("HiTEA (Germline)", "unique","HiTEA+PacB", "PacB only","PacB (Germline)"))
  #partition = res$TE
  
  ht_list = Heatmap(partition, col = structure(2:7, names = paste0(levels(as.factor(partition)) )), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    EnrichedHeatmap(m.unspl.hi,axis_name = c("","Insertion",""), col = col, name = "#reads(c)",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))),
                    column_title = "Repeat anchored \nmates",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18)) + 
    EnrichedHeatmap(m.spl.hi,axis_name = c("","Insertion",""), col = col, name = "#reads(d)",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))), 
                    column_title = "Repeat mapping \nclipped reads",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18))
  
  png(paste0(FIGDIR,"CallsCoverageComparison_GM12878_50xsc.png"),width = 400,height = 800)
  draw(ht_list, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
}

## unique calls (HiTEA based discordant reads in WGS)
if(F){
  load(paste0(DATADIR,"CoveragePlotMatrices_50Xsc.RData"))
  
  res <- as.data.frame(res)
  res$call <- gsub("Germline","HiTEA (Germline)",res$call)
  res$call <- gsub("PacB HiTEA [(]Germline[)]","PacB (Germline)",res$call)
  table(res$call)
  res <- res[res$call=="unique",]
  
  m.spl.wgs <- subset(m.spl.wgs,rownames(m.spl.wgs)%in%res$id)
  m.unspl.wgs <- subset(m.unspl.wgs,rownames(m.unspl.wgs)%in%res$id)
  
  m.spl.wgs <- makematrix(m.spl.wgs)
  m.unspl.wgs <- makematrix(m.unspl.wgs)
  maxval <- c(max(m.spl.wgs),max(m.unspl.wgs))
  maxval <- 6
  col = colorRamp2(c(0, max(maxval)), c("white","black"))
  
  partition = rep(1,nrow(m.spl.wgs))
  ht_list = Heatmap(partition, col = structure(2:7, names = paste0(levels(as.factor(partition)) )), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    EnrichedHeatmap(m.unspl.wgs,axis_name = c("","Insertion",""), col = col, name = "#reads(c)",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))),
                    column_title = "Repeat anchored \nmates",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18)) + 
    EnrichedHeatmap(m.spl.wgs,axis_name = c("","Insertion",""), col = col, name = "#reads(d)",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))), 
                    column_title = "Repeat mapping \nclipped reads",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18))
  
  png(paste0(FIGDIR,"CallsCoverageComparison_GM12878_50xsc_HiTEAuq_HiTEABased.png"),width = 400,height = 400)
  draw(ht_list, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
}
## 
## unique calls (MELT based discordant reads)
if(F){
  load(paste0(DATADIR,"CoveragePlotMatrices_50Xsc_discWGS.RData")) 
  m.spl.wgs <- subset(m.spl.wgs,rownames(m.spl.wgs)%in%res[res$call=="unique",]$id)
  m.unspl.wgs <- subset(m.unspl.wgs,rownames(m.unspl.wgs)%in%res[res$call=="unique",]$id)
  m.spl.wgs <- makematrix(m.spl.wgs)
  m.unspl.wgs <- makematrix(m.unspl.wgs)
  maxval <- c(max(m.spl.wgs),max(m.unspl.wgs))
  col = colorRamp2(c(0, max(maxval)), c("white","gray8"))
  
  res <- as.data.frame(res)
  res$call <- gsub("Germline","HiTEA (Germline)",res$call)
  res$call <- gsub("PacB HiTEA [(]Germline[)]","PacB (Germline)",res$call)
  table(res$call)
  
  partition = rep(1,nrow(m.spl.wgs))
  ht_list = Heatmap(partition, col = structure(2:7, names = paste0(levels(as.factor(partition)) )), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    EnrichedHeatmap(m.unspl.wgs,axis_name = c("","Insertion",""), col = col, name = "#reads(c)",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))),
                    column_title = "Repeat anchored \nmates",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18)) + 
    EnrichedHeatmap(m.spl.wgs,axis_name = c("","Insertion",""), col = col, name = "#reads(d)",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))), 
                    column_title = "Repeat mapping \nclipped reads",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18))
  
  png(paste0(FIGDIR,"CallsCoverageComparison_GM12878_50xsc_HiTEAuq.png"),width = 400,height = 400)
  draw(ht_list, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
}

########################################
####################  20x  ####################
if(F){
  library(EnrichedHeatmap)
  library(circlize)
  
  load(paste0(DATADIR,"CoveragePlotMatrices_20Xsc.RData"))
  table(res$call)
  m.spl.hi <- makematrix(m.spl.hi)
  m.unspl.hi <- makematrix(m.unspl.hi)
  maxval <- c(max(m.spl.hi),max(m.unspl.hi))
  col = colorRamp2(c(0, max(maxval)), c("white","black"))
  
  res <- as.data.frame(res)
  res$call <- gsub("Germline","HiTEA (Germline)",res$call)
  res$call <- gsub("PacB HiTEA [(]Germline[)]","PacB (Germline)",res$call)
  table(res$call)
  
  partition = res$call
  partition <- factor(partition,levels=c("HiTEA (Germline)", "unique","HiTEA+PacB", "PacB only","PacB (Germline)"))
  #partition = res$TE
  
  ht_list = Heatmap(partition, col = structure(2:7, names = paste0(levels(as.factor(partition)) )), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    EnrichedHeatmap(m.unspl.hi,axis_name = c("","Insertion",""), col = col, name = "#reads(c)",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))),
                    column_title = "Repeat anchored \nmates",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18)) + 
    EnrichedHeatmap(m.spl.hi,axis_name = c("","Insertion",""), col = col, name = "#reads(d)",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:7),yaxis_facing = "left",yaxis_gp = gpar(fontsize=12))), 
                    column_title = "Repeat mapping \nclipped reads",column_title_gp = gpar(fontsize = 20),axis_name_gp = gpar(fontsize = 18))
  
  png(paste0(FIGDIR,"CallsCoverageComparison_GM12878_20xsc.png"),width = 400,height = 800)
  draw(ht_list, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
}


########################################
########################################