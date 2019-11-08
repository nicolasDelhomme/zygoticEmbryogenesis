#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))

#' * Helper files
suppressMessages(source(here("UPSCb-common/src/R/plotMA.R")))
# TODO use here() everywhere

suppressMessages(source(here("UPSCb-common/src/R/volcanoPlot.R")))

#' * Graphics
pal=brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

# TODO need to use biological.R with only ZE
#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds,vst,gene_id){
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    return(
        # TODO - adjust the x, col and group (to your metadata)
        
        ggplot(bind_cols(as.data.frame(colData(dds)),
                         melt(vst[sel,])),
               aes(x=Time,y=value,col=Tissue,group=Tissue)) +
            geom_point() + geom_smooth() +
            scale_y_continuous(name="VST expression") + 
            ggtitle(label=paste("Expression for: ",gene_id))
    )
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds)){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj)
    
    if(verbose){
        message(sprintf("There are %s genes that are DE",sum(sel)))
    }
            
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        write.csv(rownames(res[sel,])[res[sel,"log2FoldChange"]>0],file.path(default_dir,paste0(default_prefix,"upgenes.csv")))
        write.csv(rownames(res[sel,])[res[sel,"log2FoldChange"]<0],file.path(default_dir,paste0(default_prefix,"downgenes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel]
        )
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel,])[res[sel,"log2FoldChange"]>0],
                dn=rownames(res[sel,])[res[sel,"log2FoldChange"]<0]))
}

#' # Analysis
#' * Data
load(here("analysis/salmon/ZE-ZF-Dataset-dds.rda"))
dds_DE_ZE_Cache <- dds
dds_DE_ZE_Cache$NGI.ID
dds_DE_ZE_Cache <- dds_DE_ZE_Cache[,!(dds_DE_ZE_Cache$NGI.ID == "P11562_148")]


#' ## Normalisation for visualisation
#' the normalisation is aware to take advantage of the model to determine the dispersion
vsd <- varianceStabilizingTransformation(dds_DE_ZE_Cache,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
vst_DE_ZE_Cache <- vst

#' ## Gene of interest
#' * 245379 
#' The gene is medium expressed and show different patterns
#' in ECM and FLM
line_plot(dds_DE_ZE_Cache,vst,"MA_99998g0010.1")

#' ## Differential Expression
dds_DE_ZE_Cache <- DESeq(dds_DE_ZE_Cache)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds_DE_ZE_Cache)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `ECM` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds_DE_ZE_Cache)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### FLM _vs._ ECM at T3


B2vsB1FMG <- extract_results(dds_DE_ZE_Cache,vst,c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                             export = FALSE,plot = FALSE,
                         default_prefix="Tissue_B2-FMG_B1-FMG__",
                         labels=paste0(colData(dds_DE_ZE_Cache)$Tissue,
                                       colData(dds_DE_ZE_Cache)$Time),
                         sample_sel=colData(dds_DE_ZE_Cache)$Tissue!="ZE")



####original, gives total bar plot, no distinction between up or down regulated.
ndeg <- lapply(2:3,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>2){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1 
    }
    return(list(length(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = FALSE,
                    export = FALSE,plot = FALSE)),vec))
})


barplot(unlist(ndeg),beside = TRUE,las=2)
names(ndeg) <- resultsNames(dds_DE_ZE_Cache)[2:10]

ndeg3_export_names <- resultsNames(dds_DE_ZE_Cache)[2:10]
ndeg3_export_names <- str_replace(ndeg3_export_names, "Time", "FMG_Time")
ndeg3_export_names <- str_replace(ndeg3_export_names, "vs_B1", "vs_B(-1)")

####Time vs B1 FMG Tissue
###differential expression based on Time only
ndeg3 <- sapply(2:10,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>2){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg3_export_names[i-1]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[2:10]
barnames <- str_replace_all(barnames,"_","")
barnames <- str_replace(barnames,"Time"," ")
barnames <- str_replace(barnames,"vs"," vs ")
barnames <- str_replace(barnames,"vs B1","vs B")
barnames <- str_c(barnames[1:9],1:9)


#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg3) <- barnames
rownames(ndeg3) <- c("all","up","dn")
barplot(ndeg3,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg3trans <- t(ndeg3)
rownames(ndeg3trans) <- barnames
colnames(ndeg3trans) <- c("all","up","dn")
barplot(ndeg3trans,beside = TRUE, las=2, horiz=F)

resultsNames(dds_DE_ZE_Cache)
B6vsB1 <- sapply(extract_results(dds_DE_ZE_Cache,vst,c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                          export = FALSE,plot = FALSE),length)
B6vsB1 <- data.frame(B6vsB1)

B7vsB1 <- sapply(extract_results(dds_DE_ZE_Cache,vst,c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                          export = FALSE,plot = FALSE),length)
#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(c("B6 vs B1"=B6vsB1[1],
                            "B7 vs B1"=B7vsB1[1]),
                       NULL,
                       fill=pal[1:2]))




ndeg4_export_names <- resultsNames(dds_DE_ZE_Cache)[12:20]
ndeg4_export_names <- str_replace(ndeg4_export_names, "Time", "ZE_Time")
ndeg4_export_names <- str_replace(ndeg4_export_names, "vs_B1", "vs_B(-1)")


####Time vs ZETissue (not adjusted for FMG vs ZE)
ndeg4 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg4_export_names[i-11]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"ZE","vs B")
barnames <- str_c(barnames[1:9],1:9)

#Differential expression ZE vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg4) <- barnames
rownames(ndeg4) <- c("all","up","dn")
barplot(ndeg4,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg4trans <- t(ndeg4)
rownames(ndeg4trans) <- barnames
colnames(ndeg4trans) <- c("all","up","dn")
barplot(ndeg4trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG






####Time vs ZETIssue (maybe adjusted for FMG vs ZE?)
##0.5 for every ZE_Time sample, 1 for ZEvsFMG Tissue
ndeg5 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-0.5
    vec[13]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg5) <- resultsNames(dds_DE_ZE_Cache)[12:20]
rownames(ndeg5) <- c("all","up","dn")
barplot(ndeg5,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg5trans <- t(ndeg5)
rownames(ndeg5trans) <- resultsNames(dds_DE_ZE_Cache)[12:20]
colnames(ndeg5trans) <- c("all","up","dn")
barplot(ndeg5trans,beside = TRUE, las=2, horiz=F)
    ###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG



##0.5 for ZEvsFMG Tissue
ndeg6 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[11]<-0.5
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg6) <- resultsNames(dds_DE_ZE_Cache)[12:20]
rownames(ndeg6) <- c("all","up","dn")
barplot(ndeg6,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg6trans <- t(ndeg6)
rownames(ndeg6trans) <- resultsNames(dds_DE_ZE_Cache)[12:20]
colnames(ndeg6trans) <- c("all","up","dn")
barplot(ndeg6trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG





ndeg7_export_names <- resultsNames(dds_DE_ZE_Cache)[12:20]
ndeg7_export_names <- str_replace(ndeg7_export_names, "Time", "ZE_vs_FMG")
ndeg7_export_names <- str_replace(ndeg7_export_names, "vs_B1", "vs_B(-1)")

###comparing ZE and FMG time points, with -1 on the previous time point of FMG
ndeg7 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[i-10]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)

        vec[i-11] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg7_export_names[i-11]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"ZE","vs B")
barnames <- str_c(barnames[1:9],1:9)

#Differential expression ZE vs FMG (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...) - looking at the tissue effect itself
#the tissue differences between FMG and ZE
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg7) <- barnames
rownames(ndeg7) <- c("all","up","dn")
barplot(ndeg7,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg7trans <- t(ndeg7)
rownames(ndeg7trans) <- barnames
colnames(ndeg7trans) <- c("all","up","dn")
barplot(ndeg7trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG







##1 for ZE and FMG Time points, -1 to previous FMG and ZE time point
ndeg8 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[i-10]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-11] <- -1     
        vec[i-1] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"ZE","vs B")
barnames <- str_c(barnames[1:9],1:9)

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg8) <- barnames
rownames(ndeg8) <- c("all","up","dn")
barplot(ndeg8,beside = TRUE, las=2, horiz=F)







#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg8trans <- t(ndeg8)
rownames(ndeg8trans) <- barnames
colnames(ndeg8trans) <- c("all","up","dn")
barplot(ndeg8trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG




#Tissue Expression Specificity using "expressionSpecificityUtility.R"
#Looking at how genes are expressed - only in one tissue? Ubiquitously? 
#make new vector, for each time point have Tissue combined with Time (eg ZE_B1, FMG_B4, etc..)
#can run simple or complete - will be running complete
TT_mat <- str_c(dds_DE_ZE_Cache$Tissue,dds_DE_ZE_Cache$Time)
colnames(vst)


source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
TT_exp.specif <- expressionSpecificity(exp.mat = vst_DE_ZE_Cache, tissue = TT_mat, mode = c("local"), output = c("complete"))
    #exp.mat is vst_DE_ZE_Cache
#tissue will be combination of time and tissue

str_subset(colnames(TT_exp.specif), "n", negate = TRUE)
colnames(TT_exp.specif)


TT_exp.specif_clipped <- TT_exp.specif[,2:21]
#now that I have the tissue specific expression, I should subset the genes - look for genes that are DE
# extract genes that are DE, and pull those genes out and make a heatmap to represent the genes
pc.TT_exp.specif <- prcomp(t(vst_DE_ZE_Cache))
got_pc <- get_pca(pc.TT_exp.specif, element = c("var", "ind"))

genes <- sort(got_pc$contrib[,1],decreasing = TRUE)
genes_PC1_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1]]
genes <- sort(got_pc$contrib[,2],decreasing = TRUE)
genes_PC2_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,2],decreasing = TRUE))>=50)[1]]
genes <- sort(got_pc$contrib[,3],decreasing = TRUE)
genes_PC3_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,3],decreasing = TRUE))>=50)[1]]

length(genes_PC1_50) #1239
length(genes_PC2_50) #1253
length(genes_PC3_50) #1358

genes_PC1_50_names <- names(genes_PC1_50)
genes_PC2_50_names <- names(genes_PC2_50)
genes_PC3_50_names <- names(genes_PC3_50)

#background <- rownames(vst_DE_ZE_Cache)
#background <- str_replace(background,"[.]1","")

#genes_PC1_50_names <- str_replace(genes_PC1_50_names,"[.]1","")
#genes_PC2_50_names <- str_replace(genes_PC2_50_names,"[.]1","")
#genes_PC3_50_names <- str_replace(genes_PC3_50_names,"[.]1","")



###fix these heatmaps, need to normalize them?

#PC1 Top 50%
TT_exp.specif.PC1 <- TT_exp.specif_clipped[genes_PC1_50_names,]
heatmap.2(TT_exp.specif.PC1,
          col = hpal)

#PC2 Top 50%
TT_exp.specif.PC2 <- TT_exp.specif_clipped[genes_PC2_50_names,]
heatmap.2(TT_exp.specif.PC2,
          col = hpal)

#PC3 Top 50%
TT_exp.specif.PC3 <- TT_exp.specif_clipped[genes_PC3_50_names,]
heatmap.2(TT_exp.specif.PC3,
          col = hpal)

#DE Expressed Genes FMG vs Time
#DE Expressed Genes ZE vs Time
#DE Expressed Genes ZE vs FMG






#relevel
dds_DE_ZE_Cache_Releveled <- dds_DE_ZE_Cache
relevel(dds_DE_ZE_Cache_Releveled$Tissue, "ZE")
dds_DE_ZE_Cache_Releveled$Tissue <- relevel(dds_DE_ZE_Cache_Releveled$Tissue, "ZE")

dds_DE_ZE_Cache_Releveled <- DESeq(dds_DE_ZE_Cache_Releveled)
resultsNames(dds_DE_ZE_Cache_Releveled)



#NDEG Releveled FMG vs Time
{
ndeg_relevel_FMG_Time <- sapply(2:10,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i]<-1
    if(i>2){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-1]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[2:10]
barnames <- str_replace_all(barnames,"_","")
barnames <- str_replace(barnames,"Time"," ")
barnames <- str_replace(barnames,"vs"," vs ")
barnames <- str_replace(barnames,"vs B1","vs B")
barnames <- str_c(barnames[1:9],1:9)


#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_relevel_FMG_Time) <- barnames
rownames(ndeg_relevel_FMG_Time) <- c("all","up","dn")
barplot(ndeg_relevel_FMG_Time,beside = TRUE, las=2, horiz=F)
}


#NDEG Releveled ZE vs Time
{
ndeg_relevel_ZE_Time <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg4_export_names[i-11]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"ZE","vs B")
barnames <- str_c(barnames[1:9],1:9)

#Differential expression ZE vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_relevel_ZE_Time) <- barnames
rownames(ndeg_relevel_ZE_Time) <- c("all","up","dn")
barplot(ndeg_relevel_ZE_Time,beside = TRUE, las=2, horiz=F)
}

#NDEG Releveled ZE vs FMG
{
ndeg_relevel_ZE_vs_FMG <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i]<-1
    vec[i-10]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        
        vec[i-11] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg7_export_names[i-11]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"FMG","vs B")
barnames <- str_c(barnames[1:9],1:9)

#Differential expression ZE vs FMG (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...) - looking at the tissue effect itself
#the tissue differences between FMG and ZE
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_relevel_ZE_vs_FMG) <- barnames
rownames(ndeg_relevel_ZE_vs_FMG) <- c("all","up","dn")
barplot(ndeg_relevel_ZE_vs_FMG,beside = TRUE, las=2, horiz=F)
}

#NDEG Releveled ZE vs FMG v2 #REAL TISSUE EFFECT MAYBE?
{
    ndeg_relevel_ZE_vs_FMG.v2 <- sapply(12:20,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
        vec[i]<-1
        vec[i-10]<-1
        if(i>12){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            
            vec[i-11] <- -1
            vec[i-1] <- -1
        }
        sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst,vec,verbose = TRUE,
                               export = FALSE,plot = FALSE, default_prefix = ndeg7_export_names[i-11]),length)
    })
    
    barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[12:20]
    barnames <- str_replace(barnames,"Time","")
    barnames <- str_replace(barnames,".Tissue"," ")
    barnames <- str_replace(barnames,"FMG","vs B")
    barnames <- str_c(barnames[1:9],1:9)
    
    #Differential expression ZE vs FMG (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...) - looking at the tissue effect itself
    #the tissue differences between FMG and ZE
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_relevel_ZE_vs_FMG.v2) <- barnames
    rownames(ndeg_relevel_ZE_vs_FMG.v2) <- c("all","up","dn")
    barplot(ndeg_relevel_ZE_vs_FMG.v2,beside = TRUE, las=2, horiz=F)
}





#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("B2-B3"=ZEB2B3$all,
                            "B4-B6"=ZEB4B5B6$all,
                            "B7-B10"=ZEB7B8B910$all),
                       NULL,
                       fill=pal[1:3]))




#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


