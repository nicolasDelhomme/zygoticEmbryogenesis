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
                         data.frame(value=vst[sel,])),
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
#' Caching the data so that it won't be overwritten by other scripts.
load(here("analysis/salmon/ZE-ZF-Dataset-dds.rda"))

dds_DE_ZE_Cache <- dds.ze
#dds_DE_ZE_Cache <- dds

dds_DE_ZE_Cache$NGI.ID
dds_DE_ZE_Cache <- dds_DE_ZE_Cache[,!(dds_DE_ZE_Cache$NGI.ID == "P11562_148")]


#' ## Normalisation for visualisation
#' the normalisation is aware to take advantage of the model to determine the dispersion
vsd <- varianceStabilizingTransformation(dds_DE_ZE_Cache,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
vst_DE_ZE_Cache <- vst

#' ## Lineplots
rownames(dds_DE_ZE_Cache)

line_plot(dds_DE_ZE_Cache,vst,"MA_99998g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_100150g0010.1")

line_plot(dds_DE_ZE_Cache,vst,"MA_10101060g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10106144g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10098066g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10098872g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10107300g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10099160g0010.1")
line_plot(dds_DE_ZE_Cache,vst,"MA_10086594g0010.1")


dds_DE_ZE_Cache$NGI.ID
#' ## Differential Expression
dds_DE_ZE_Cache <- DESeq(dds_DE_ZE_Cache)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds_DE_ZE_Cache)

#' The model used is:
#' 
#' Tissue * Time
resultsNames(dds_DE_ZE_Cache)

#relevel
dds_DE_ZE_Cache_Releveled <- dds_DE_ZE_Cache
relevel(dds_DE_ZE_Cache_Releveled$Tissue, "ZE")
dds_DE_ZE_Cache_Releveled$Tissue <- relevel(dds_DE_ZE_Cache_Releveled$Tissue, "ZE")

dds_DE_ZE_Cache_Releveled <- DESeq(dds_DE_ZE_Cache_Releveled)
resultsNames(dds_DE_ZE_Cache_Releveled)




#' ## Results

# Time effect with the reference level of FMG
ndeg_Time_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
barnames <- str_replace_all(barnames,"_","")
barnames <- str_replace(barnames,"Time"," ")
barnames <- str_replace(barnames,"vs"," vs ")
barnames <- str_replace(barnames,"vs B1","vs B")
barnames <- str_c(barnames[1:9],1:9)
}
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_Time_ZF) <- barnames
rownames(ndeg_Time_ZF) <- c("all","up","dn")

barplot(ndeg_Time_ZF,beside = TRUE, las=2, horiz=F)
}

# Interaction Term
ndeg_TissueInteraction_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
    barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
    barnames <- str_replace_all(barnames,"_","")
    barnames <- str_replace(barnames,"Time"," ")
    barnames <- str_replace(barnames,"vs"," vs ")
    barnames <- str_replace(barnames,"vs B1","vs B")
    barnames <- str_c(barnames[1:9],1:9)
}    
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_TissueInteraction_ZF) <- barnames
rownames(ndeg_TissueInteraction_ZF) <- c("all","up","dn")  

barplot(ndeg_TissueInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

#Main Effect + Interaction Term
#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_MainInteraction_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    vec[i] <- 1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
    barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
    barnames <- str_replace_all(barnames,"_","")
    barnames <- str_replace(barnames,"Time"," ")
    barnames <- str_replace(barnames,"vs"," vs ")
    barnames <- str_replace(barnames,"vs B1","vs B")
    barnames <- str_c(barnames[1:9],1:9)
    
}    
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_MainInteraction_ZF) <- barnames
rownames(ndeg_MainInteraction_ZF) <- c("all","up","dn")

barplot(ndeg_MainInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

# Interaction Term ALONE
ndeg_TissueInteractionAlone_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1

    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TissueInteractionAlone_ZF) <- barnames
    rownames(ndeg_TissueInteractionAlone_ZF) <- c("all","up","dn")  
    
    barplot(ndeg_TissueInteractionAlone_ZF,beside = TRUE, las=2, horiz=F)
}

#Main Effect + Interaction Term ALONE
#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_MainInteractionAlone_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    vec[i] <- 1

    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_MainInteractionAlone_ZF) <- barnames
    rownames(ndeg_MainInteractionAlone_ZF) <- c("all","up","dn")
    
    barplot(ndeg_MainInteractionAlone_ZF,beside = TRUE, las=2, horiz=F)
}

#DEVEL BRANCH
#Main + Timepoints, -1 to previous
ndeg_MainTime_ZF <- sapply(3:11,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        vec[2] <- 1
        vec[i] <- 1
        if(i>3){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            vec[i-1] <- -1        
        }
        sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                               export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
    })
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_MainTime_ZF) <- barnames
    rownames(ndeg_MainTime_ZF) <- c("all","up","dn")
    
    barplot(ndeg_MainTime_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints - Main, -1 to previous
ndeg_TimenoMain_ZF <- sapply(3:11,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        vec[2] <- -1
        vec[i] <- 1
        if(i>3){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            vec[i-1] <- -1        
        }
        sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                               export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
    })
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimenoMain_ZF) <- barnames
    rownames(ndeg_TimenoMain_ZF) <- c("all","up","dn")
    
    barplot(ndeg_TimenoMain_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction, -1 to previous on Time
ndeg_TimeInteraction_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    vec[i+9] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeInteraction_ZF) <- barnames
    rownames(ndeg_TimeInteraction_ZF) <- c("all","up","dn")
    
    barplot(ndeg_TimeInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction, -1 to previous on Time AND Interaction
ndeg_TimeInteractionDoubleContrast_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    vec[i+9] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1
        vec[i+8] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeInteractionDoubleContrast_ZF) <- barnames
    rownames(ndeg_TimeInteractionDoubleContrast_ZF) <- c("all","up","dn")
    
    barplot(ndeg_TimeInteractionDoubleContrast_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction + Main Tissue Effect, -1 to previous on Time
ndeg_TimeInteractionMain_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    vec[i+9] <- 1
    vec[2] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeInteractionMain_ZF) <- barnames
    rownames(ndeg_TimeInteractionMain_ZF) <- c("all","up","dn")
    
    barplot(ndeg_TimeInteractionMain_ZF,beside = TRUE, las=2, horiz=F)
}


#Now do a re-level
resultsNames(dds_DE_ZE_Cache_Releveled)
#' Time effect with the reference level of FMG
ndeg_releveled_Time_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i]<-1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
    barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
    barnames <- str_replace_all(barnames,"_","")
    barnames <- str_replace(barnames,"Time"," ")
    barnames <- str_replace(barnames,"vs"," vs ")
    barnames <- str_replace(barnames,"vs B1","vs B")
    barnames <- str_c(barnames[1:9],1:9)
}
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_releveled_Time_ZF) <- barnames
rownames(ndeg_releveled_Time_ZF) <- c("all","up","dn")

barplot(ndeg_releveled_Time_ZF,beside = TRUE, las=2, horiz=F)
}

# Interaction Term
ndeg_releveled_TissueInteraction_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
    barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
    barnames <- str_replace_all(barnames,"_","")
    barnames <- str_replace(barnames,"Time"," ")
    barnames <- str_replace(barnames,"vs"," vs ")
    barnames <- str_replace(barnames,"vs B1","vs B")
    barnames <- str_c(barnames[1:9],1:9)
}    
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_releveled_TissueInteraction_ZF) <- barnames
rownames(ndeg_releveled_TissueInteraction_ZF) <- c("all","up","dn")  

barplot(ndeg_releveled_TissueInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

#Main Effect + Interaction Term
#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_releveled_MainInteraction_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[2] <- 1
    vec[i] <- 1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
#Barplot Preprocessing
{
    barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
    barnames <- str_replace_all(barnames,"_","")
    barnames <- str_replace(barnames,"Time"," ")
    barnames <- str_replace(barnames,"vs"," vs ")
    barnames <- str_replace(barnames,"vs B1","vs B")
    barnames <- str_c(barnames[1:9],1:9)
    
}    
#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg_releveled_MainInteraction_ZF) <- barnames
rownames(ndeg_releveled_MainInteraction_ZF) <- c("all","up","dn")

barplot(ndeg_releveled_MainInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

# Interaction Term ALONE
ndeg_releveled_TissueInteractionAlone_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
#    if(i>12){
#        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
#        vec[i-10] <- -1
#    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_TissueInteractionAlone_ZF) <- barnames
    rownames(ndeg_releveled_TissueInteractionAlone_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_TissueInteractionAlone_ZF,beside = TRUE, las=2, horiz=F)
}

#Main Effect + Interaction Term ALONE
#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_releveled_MainInteractionAlone_ZF <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[2] <- 1
    vec[i] <- 1

    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_MainInteractionAlone_ZF) <- barnames
    rownames(ndeg_releveled_MainInteractionAlone_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_MainInteractionAlone_ZF,beside = TRUE, las=2, horiz=F)
}


#DEVEL BRANCH
#Main + Timepoints, -1 to previous
ndeg_releveled_MainTime_ZF <- sapply(3:11,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        vec[2] <- 1
        vec[i] <- 1
        if(i>3){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            vec[i-1] <- -1        
        }
        sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                               export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
    })
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_MainTime_ZF) <- barnames
    rownames(ndeg_releveled_MainTime_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_MainTime_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints - Main, -1 to previous
ndeg_releveled_TimenoMain_ZF <- sapply(3:11,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
        vec[2] <- -1
        vec[i] <- 1
        if(i>3){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            vec[i-1] <- -1        
        }
        sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                               export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
    })
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_TimenoMain_ZF) <- barnames
    rownames(ndeg_releveled_TimenoMain_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_TimenoMain_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction, -1 to previous on Time
ndeg_releveled_TimeInteraction_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
    vec[i+9] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_TimeInteraction_ZF) <- barnames
    rownames(ndeg_releveled_TimeInteraction_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_TimeInteraction_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction, -1 to previous on Time AND Interaction
ndeg_releveled_TimeInteractionDoubleContrast_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
    vec[i+9] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1
        vec[i+8] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_TimeInteractionDoubleContrast_ZF) <- barnames
    rownames(ndeg_releveled_TimeInteractionDoubleContrast_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_TimeInteractionDoubleContrast_ZF,beside = TRUE, las=2, horiz=F)
}

#Timepoints + Interaction + Main Tissue Effect, -1 to previous on Time
ndeg_releveled_TimeInteractionMain_ZF <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
    vec[i+9] <- 1
    vec[2] <- 1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_releveled_TimeInteractionMain_ZF) <- barnames
    rownames(ndeg_releveled_TimeInteractionMain_ZF) <- c("all","up","dn")
    
    barplot(ndeg_releveled_TimeInteractionMain_ZF,beside = TRUE, las=2, horiz=F)
}

#final ndeg values
ndeg_Time_ZF
ndeg_TissueInteraction_ZF
ndeg_MainInteraction_ZF
ndeg_MainTime_ZF
ndeg_TimenoMain_ZF
ndeg_TimeInteraction_ZF
ndeg_TimeInteractionDoubleContrast_ZF

ndeg_TimeInteractionMain_ZF

ndeg_TissueInteractionAlone_ZF
ndeg_MainInteractionAlone_ZF

#final ndeg values releveled
ndeg_releveled_Time_ZF
ndeg_releveled_TissueInteraction_ZF
ndeg_releveled_MainInteraction_ZF
ndeg_releveled_MainTime_ZF
ndeg_releveled_TimenoMain_ZF
ndeg_releveled_TimeInteraction_ZF
ndeg_releveled_TimeInteractionDoubleContrast_ZF

ndeg_releveled_TimeInteractionMain_ZF

ndeg_releveled_TissueInteractionAlone_ZF
ndeg_releveled_MainInteractionAlone_ZF



# Time effect with the reference level of FMG
ndeg_Test <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <-1
    vec[2] <- 1
    vec[i-9] <- 1
    if(i>13){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1
        vec[i-10] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
    }
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_Time_ZF) <- barnames
    rownames(ndeg_Time_ZF) <- c("all","up","dn")
    
    barplot(ndeg_Time_ZF,beside = TRUE, las=2, horiz=F)
}

# Time effect with the reference level of FMG
ndeg_Test <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    vec[2] <- 1
    # vec[i-9] <- -1 IS REMOVING THE EFFECT OF FMG, THEREFORE WE SEE THE EXPRESSION OF ZE TISSUE AT THAT TIME
    #WE AREN'T JUST REMOVING THE TIME EFFECT, WE ARE REMOVING THE TISSUE EFFECT RELATED WITH THAT TIME INFORMATION
    vec[i-9] <- -1
    if(i>13){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1
        vec[i-10] <- -1
    }
    print(vec)
    print(resultsNames(dds_DE_ZE_Cache))
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
    }
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_Time_ZF) <- barnames
    rownames(ndeg_Time_ZF) <- c("all","up","dn")
    
    barplot(ndeg_Time_ZF,beside = TRUE, las=2, horiz=F)
}

#This gives us B2vsB1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,c("Time","B1","B2"),verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 3048 DE genes

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[3] <- 1 #B2vsB1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 3048 DE genes

resultsNames(dds_DE_ZE_Cache)
#THIS CONFIRMS THAT WE ARE DOING THE TIME CONTRAST CORRECTLY (-1 for PREVIOUS)
#This gives us B4vsB3, in FMG only
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,c("Time","B3","B4"),verbose = TRUE,
                export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
#This gives me 724 correct in NDEG normal time

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[5] <- 1 #B4
test_vec[4] <- -1 #B3
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
#This gives me 724 - correct in NDEG normal time

#effectively, B4_vs_B3


#following along to "The effect of treatment in mutant"
test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[3] <- 1 #B2vsB1
test_vec[12] <- 1 #TissueEffect
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 3050 DE genes

#something wrong, 3048 DE in FMG, whereas there is 3050 DE in ZE (effect + interaction)????






test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[5] <- 1 #B4
test_vec[4] <- -1 #B3
test_vec[14] <- 1 #int B4
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 2022 DE genes - matches "Tissue Interaction" plot, this should be the effect in ZE

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[6] <- 1 #B5
test_vec[5] <- -1 #B4
test_vec[15] <- 1 #int B5
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 2377 DE genes

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[7] <- 1 #B6
test_vec[6] <- -1 #B5
test_vec[16] <- 1 #int B6
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 2916 DE genes

#looking at above, without interaction
test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[7] <- 1 #B6
test_vec[6] <- -1 #B5
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 816 DE genes - same as in NDEG normal time plot

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[16] <- 1 #int B6
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
#just interaction from above DE genes are 1595

#difference between time points is:
# 1) ZE_vs_FMG
# 2) above and TissueZE.Time(#)

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[2] <- 1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 0 DE - good to see

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[2] <- 1
test_vec[12] <- 1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 0 DE - very good to see

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[2] <- 1
test_vec[13] <- 1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 0 DE - very nice

test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[2] <- 1
test_vec[14] <- 1
TIMETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                            export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])
# 6064 DE


#### Time vs FMG (baseline) - put that here
ndeg_TimeDifferenceFMG <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    if(i>3){
        vec[i-1] <- -1
    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeDifferenceFMG) <- barnames
    rownames(ndeg_TimeDifferenceFMG) <- c("all","up","dn")
    
    barplot(ndeg_TimeDifferenceFMG,beside = TRUE, las=2, horiz=F, ylim = c(0,8000))
}


#### Time vs ZE (baseline + interaction)
#B2_vs_B1 = 1, GenoB2 = 1
#B3_vs_B1 = 1, B2_vs_B1 = -1, GenoB3 = 1
ndeg_TimeDifferenceZE <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i] <- 1
    vec[i+9] <- 1
    if(i>3){
        vec[i-1] <- -1
    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeDifferenceZE) <- barnames
    rownames(ndeg_TimeDifferenceZE) <- c("all","up","dn")
    
    barplot(ndeg_TimeDifferenceZE,beside = TRUE, las=2, horiz=F, ylim = c(0,8000))
}

ndeg_TimeDifferenceZEReleveled <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
    if(i>3){
        vec[i-1] <- -1
    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TimeDifferenceZEReleveled) <- barnames
    rownames(ndeg_TimeDifferenceZEReleveled) <- c("all","up","dn")
    
    barplot(ndeg_TimeDifferenceZEReleveled,beside = TRUE, las=2, horiz=F, ylim = c(0,8000))
}





#### now make the proper thing
# Difference Between Tissue Types - Just Interactions <- this is the same as MainInteractionAlone plot!
ndeg_TissueDifference <- sapply(2:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    if(i>2){
        vec[i+9] <- 1
    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[2:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:10],1:10)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TissueDifference) <- barnames
    rownames(ndeg_TissueDifference) <- c("all","up","dn")
    
    barplot(ndeg_TissueDifference,beside = TRUE, las=2, horiz=F, ylim = c(0,12000))
}
#### This is the DE gene expression based on tissue


# Difference Between Tissue Types - Just Interactions
ndeg_TissueDifferenceMark2 <- sapply(2:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    if(i==2){
        vec[2] <- 1
    }else if(i==3){
        vec[2] <- 1
        vec[i+9] <- 1
    }else if(i>3){
        vec[i+8] <- 1
        vec[i+9] <- 1
    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                          export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache)[2:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:10],1:10)
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_TissueDifferenceMark2) <- barnames
    rownames(ndeg_TissueDifferenceMark2) <- c("all","up","dn")
    
    barplot(ndeg_TissueDifferenceMark2,beside = TRUE, las=2, horiz=F, ylim = c(0,12000))
}









#subset only time point B7 to see if we get this differential expression



#only looking at time point B7

dds_subset <- dds_DE_ZE_Cache

dds_subset <- dds_subset[,(dds_subset$Time == "B7")]
dds_subset$Time <- droplevels(dds_subset$Time)

colData(dds_subset)
design(dds_subset) <- ~Tissue

dds_subset <- DESeq(dds_subset)
resultsNames(dds_subset)
vst_subset <- vst_DE_ZE_Cache
vst_subset <- vst_subset[,dds_subset$NGI.ID]

DE_B7_subset <- sapply(1:1,function(i){
    vec <- rep(0,length(resultsNames(dds_subset)))
    vec[i] <- 1
    #    if(i>3){
    #        vec[i-1] <- -1
    #    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_subset,vst_subset,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

res <- results(dds_subset)

res2 <- res[abs(res$log2FoldChange) >= 0.5 & !is.na(res$padj) & res$padj <= 0.01 & res$baseMean >= 1000,]


res2[order(res2$log2FoldChange,decreasing = TRUE)[1],]
line_plot(dds_subset,vst_subset, "MA_138039g0010.1")
line_plot(dds_DE_ZE_Cache,vst_DE_ZE_Cache, "MA_138039g0010.1")


    
    
for(i in 1:10){ 
    print(rownames(res2[order(res2$log2FoldChange,decreasing = TRUE)[i],]))
    print(res2[order(res2$log2FoldChange,decreasing = TRUE)[i],]$baseMean)
    
    gname <- rownames(res2[order(res2$log2FoldChange,decreasing = TRUE)[i],])

    #p <- line_plot(dds_DE_ZE_Cache,vst_DE_ZE_Cache, gname)
    #plot(p)
    }



# nice looking plot of this gene, number 6 in loop above
# "MA_191819g0010.1"
# "MA_427046g0010.1"
# "MA_22203g0020.1"


#test other genes to make sure l2fc is correct

























dds_subset <- dds_DE_ZE_Cache

dds_subset <- dds_subset[,!(dds_subset$Time == "B1")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B2")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B3")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B4")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B5")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B8")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B9")]
dds_subset <- dds_subset[,!(dds_subset$Time == "B10")]
dds_subset$Time
dds_subset$Time <- droplevels(dds_subset$Time)

dds_subset <- DESeq(dds_subset)

resultsNames(dds_subset)
results(dds_subset)
design(dds_subset)

res.b67 <- results(dds_subset)

res2.b67 <- res[abs(res$log2FoldChange) >= 0.5 & !is.na(res$padj) & res$padj <= 0.01 & res$baseMean >= 1000,]
res2.b67

vst_subset <- vst_DE_ZE_Cache
vst_subset <- vst_subset[,dds_subset$NGI.ID]
line_plot(dds_subset,vst_subset, "MA_191819g0010.1")

res2["MA_191819g0010.1",]
#2.80 L2FC
res2["MA_427046g0010.1",]
#2.24 L2FC

res.tis <- results(dds_subset, name = "Tissue_ZE_vs_FMG")
res.tim <- results(dds_subset, name = "Time_B7_vs_B6")
res.tt <- results(dds_subset, name = "TissueZE.TimeB7")
res.gene <- results(dds_subset, contrast = c(0,1,0,1))

res2.tis <- res.tis[abs(res.tis$log2FoldChange) >= 0.5 & !is.na(res.tis$padj) & res.tis$padj <= 0.01 & res.tis$baseMean >= 1000,]
res2.tim <- res.tim[abs(res.tim$log2FoldChange) >= 0.5 & !is.na(res.tim$padj) & res.tim$padj <= 0.01 & res.tim$baseMean >= 1000,]
res2.tt <- res.tt[abs(res.tt$log2FoldChange) >= 0.5 & !is.na(res.tt$padj) & res.tt$padj <= 0.01 & res.tt$baseMean >= 1000,]

res.tis["MA_191819g0010.1",]$log2FoldChange
res.tim["MA_191819g0010.1",]$log2FoldChange
res.tt["MA_191819g0010.1",]$log2FoldChange
res.gene["MA_191819g0010.1",]
res.gene["MA_427046g0010.1",]
res.gene["MA_22203g0020.1",]
#2.98 l2fc

resultsNames(dds_subset)

#Tissue_ZE_vs_FMG
#The difference between ZE and FMG at the baseline (B6) time point
DE_difference_betweeN_tissues <- sapply(3:3,function(i){
    vec <- rep(0,length(resultsNames(dds_subset)))
    vec[i-1] <- 1
    #    if(i>3){
    #        vec[i-1] <- -1
    #    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_subset,vst_subset,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})



#the effect of time in reference level tissue (FMG)
DE_main_effect <- sapply(3:3,function(i){
    vec <- rep(0,length(resultsNames(dds_subset)))
    vec[i] <- 1
#    if(i>3){
#        vec[i-1] <- -1
#    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_subset,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

#The effect of time in ZE tissue
DE_effect_of_time_ZE <- sapply(3:3,function(i){
    vec <- rep(0,length(resultsNames(dds_subset)))
    vec[i] <- 1
    vec[i+1] <- 1
    #    if(i>3){
    #        vec[i-1] <- -1
    #    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    sapply(extract_results(dds_subset,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})



#The difference between ZE and FMG at the baseline (B7) time point
#ZE vs FMG (referenced against B6 - the main), TissueZE.B7 (the interaction term)
DE_difference_betweeN_tissues_at_time <- sapply(3:3,function(i){
    vec <- rep(0,length(resultsNames(dds_subset)))
    vec[2] <- 1
    vec[4] <- 1
    #    if(i>3){
    #        vec[i-1] <- -1
    #    }
    print(vec)
    #    if(i>12){
    #        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
    #        vec[i-10] <- -1
    #    }
    res.gene <- results(dds_subset, contrast = vec)
    print(res.gene["MA_191819g0010.1",])
#    sapply(extract_results(dds_subset,vst_DE_ZE_Cache,vec,verbose = TRUE,
#                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
DE_difference_betweeN_tissues_at_time
#B6 vs B7
#9963 DE genes


#############THIS IS RIGHT
#############THIS IS RIGHT
#############THIS IS RIGHT


##test with the full model
#Main Effect + Interaction Term ALONE
#ie what is the difference between genotypes ZE and FMG at the different time points
l2fc_gtest_1 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    vec[i] <- 1
 #   vec[i-9] <- -1
    
    print(vec)
    res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    print(res.gene["MA_191819g0010.1",])
    
    return(res.gene["MA_191819g0010.1",])
#    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
#                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

l2fc_gtest_2 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    vec[i] <- 1
    #   vec[i-9] <- -1
    
    print(vec)
    res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    print(res.gene["MA_427046g0010.1",])
    
    return(res.gene["MA_427046g0010.1",])
    #    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
    #                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

l2fc_gtest_3 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1
    vec[i] <- 1
    #   vec[i-9] <- -1
    
    print(vec)
    res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    print(res.gene["MA_22203g0020.1",])
    
    return(res.gene["MA_22203g0020.1",])
    #    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
    #                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

res.gene["MA_191819g0010.1",]
res.gene["MA_427046g0010.1",]
res.gene["MA_22203g0020.1",]

Reduce(rbind,l2fc_gtest_1)
Reduce(rbind,l2fc_gtest_2)
Reduce(rbind,l2fc_gtest_3)

#THIS ONE IS RIGHT - B7 is at 2.80 L2FC!!
#SIMPLY ADD TISSUE VS TISSUE AND THE INTERACTION TERM TO GET IT
#For the next one, at B7 it should be 2.24


#export names
ndeg_time_exnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
ndeg_time_exnames <- str_replace(ndeg_time_exnames, "Time", "FMG_Time")
ndeg_time_exnames <- str_replace(ndeg_time_exnames, "vs_B1", "vs_B")
ndeg_time_exnames <- str_c(ndeg_time_exnames, 1:9)

ndeg_time_ze_exnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
ndeg_time_ze_exnames <- str_replace(ndeg_time_ze_exnames, "Time", "ZE_Time")
ndeg_time_ze_exnames <- str_replace(ndeg_time_ze_exnames, "vs_B1", "vs_B")
ndeg_time_ze_exnames <- str_c(ndeg_time_ze_exnames, 1:9)

ndeg_tissue <- resultsNames(dds_DE_ZE_Cache)[2:11]
ndeg_tissue <- str_c(ndeg_tissue[1], "_B",1:10)

####THIS IS RIGHT
####THIS IS RIGHT
####THIS IS RIGHT

#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_main_add_interaction <- sapply(11:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    if(i == 11){
        vec[2] <- 1
    }else{
    vec[2] <- 1
    vec[i] <- 1
    #   vec[i-9] <- -1
    }
    print(vec)
    #res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    #print(res.gene["MA_191819g0010.1",])
    
    #return(res.gene["MA_191819g0010.1",]$log2FoldChange)
        sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                               export = TRUE,plot = FALSE, default_dir = here("analysis/DE/ZF_TissueEffect"), default_prefix = ndeg_tissue[i-10]),length)
})
#Barplot TISSUE ZE VS FMG - DIFFERENCE BETWEEN TISSUES
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[2:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
        barnames <- str_c("B",1:10)
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_main_add_interaction) <- barnames
    rownames(ndeg_main_add_interaction) <- c("all","up","dn")
    
    barplot(ndeg_main_add_interaction,beside = TRUE, las=2, horiz=F,ylim = c(0,12000))
}




#############THIS IS RIGHT
#############THIS IS RIGHT
#############THIS IS RIGHT

ndeg_time <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        vec[i] <- 1
        if(i>3){
        vec[i-1] <- -1
        }
     print(vec)
    #res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    #print(res.gene["MA_191819g0010.1",])
    
    #return(res.gene["MA_191819g0010.1",]$log2FoldChange)
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_dir = here("analysis/DE/ZF_FMG_Time"), default_prefix = ndeg_time_exnames[i-2]),length)
})
#Barplot Time FMG
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)
        
        
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_time) <- barnames
    rownames(ndeg_time) <- c("all","up","dn")
    
    barplot(ndeg_time,beside = TRUE, las=2, horiz=F,ylim = c(0,12000))
}


###Should probably relevel this one to get time and ze
ndeg_time_ze <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
    vec[i] <- 1
            if(i>3){
            vec[i-1] <- -1
            }
    print(vec)
    #res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    #print(res.gene["MA_191819g0010.1",])
    
    #return(res.gene["MA_191819g0010.1",]$log2FoldChange)
    sapply(extract_results(dds_DE_ZE_Cache_Releveled,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_dir = here("analysis/DE/ZF_ZE_Time"), default_prefix = ndeg_time_ze_exnames[i-2]),length)
})
#Barplot Time ZE
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[3:11]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:9],1:9)

    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_time_ze) <- barnames
    rownames(ndeg_time_ze) <- c("all","up","dn")
    
    barplot(ndeg_time_ze,beside = TRUE, las=2, horiz=F, ylim = c(0,12000))
}






DE_FMG_B3_B4 <- list.files(here("analysis/DE/ZF_FMG_Time"), 
                         recursive = TRUE, 
                         pattern = "",
                         full.names = TRUE)

DE_FMG_B3_B4 <- str_subset(DE_FMG_B3_B4, "Time_B4_vs_B3genes.csv")
DE_FMG_B3_B4 <- sapply(1:length(DE_FMG_B3_B4),function(i){
    
    GeneL <- read.csv(DE_FMG_B3_B4[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})
#FMG B3

DE_FMG_B7_B6 <- list.files(here("analysis/DE/ZF_FMG_Time"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)

DE_FMG_B7_B6 <- str_subset(DE_FMG_B7_B6, "Time_B7_vs_B6genes.csv")
DE_FMG_B7_B6 <- sapply(1:length(DE_FMG_B7_B6),function(i){
    
    GeneL <- read.csv(DE_FMG_B7_B6[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})
#FMG B7



DE_ZE_B3_B4 <- list.files(here("analysis/DE/ZF_ZE_Time"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)

DE_ZE_B3_B4 <- str_subset(DE_ZE_B3_B4, "Time_B4_vs_B3genes.csv")
DE_ZE_B3_B4 <- sapply(1:length(DE_ZE_B3_B4),function(i){
    
    GeneL <- read.csv(DE_ZE_B3_B4[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})
#ZE B3

DE_ZE_B7_B6 <- list.files(here("analysis/DE/ZF_ZE_Time"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)

DE_ZE_B7_B6 <- str_subset(DE_ZE_B7_B6, "Time_B7_vs_B6genes.csv")
DE_ZE_B7_B6 <- sapply(1:length(DE_ZE_B7_B6),function(i){
    
    GeneL <- read.csv(DE_ZE_B7_B6[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})
#ZE B7





#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("FMG_B3_vs_B4"=DE_FMG_B3_B4,
                            "FMG_B6_vs_B7"=DE_FMG_B7_B6,
                            "ZE_B3_vs_B4"=DE_ZE_B3_B4,
                            "ZE_B6_vs_B7"=DE_ZE_B7_B6),
                       
                       NULL,
                       fill=pal[1:4]))
venn_ZF_list <- list("FMG_B3_vs_B4"=DE_FMG_B3_B4,
                     "FMG_B6_vs_B7"=DE_FMG_B7_B6,
                     "ZE_B3_vs_B4"=DE_ZE_B3_B4,
                     "ZE_B6_vs_B7"=DE_ZE_B7_B6)

venn_ZF_list <- venn(venn_ZF_list, show.plot = FALSE)





ndeg_time[1,3]
ndeg_time[1,6]

ndeg_time_ze[1,3]
ndeg_time_ze[1,6]




















#B4vsB1 - B3-B1 = B4vsB3, with the addition of ZE interacting with B4 we get the Effect of Time in ZE
test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[5] <- 1 #B4
test_vec[4] <- -1 #B3
test_vec[14] <- 1 #ZE interacting with B4
#USING TEST_VEC THIS IS ACCURATE TO MAIN + INTERACTION ALONE
TISSUETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                              export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])

#WITH TIME, WHAT IS THE DIFFERENCE BETWEEN GENOTYPES
test_vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
test_vec[2] <- 1
test_vec[14] <- 1
TISSUETEST <- extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,test_vec,verbose = TRUE,
                              export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2])







#testing a deletion of B1, B2 and B3 because maybe its messing up everything down the line as the tissues split


dds_deleted_seedtissue <- dds_DE_ZE_Cache

dds_deleted_seedtissue <- dds_deleted_seedtissue[,!(dds_deleted_seedtissue$Time == "B1")]
dds_deleted_seedtissue <- dds_deleted_seedtissue[,!(dds_deleted_seedtissue$Time == "B2")]
dds_deleted_seedtissue <- dds_deleted_seedtissue[,!(dds_deleted_seedtissue$Time == "B3")]
dds_deleted_seedtissue$Time
dds_deleted_seedtissue$Time <- droplevels(dds_deleted_seedtissue$Time)
dds_deleted_seedtissue$Time

dds_deleted_seedtissue <- DESeq(dds_deleted_seedtissue)
vst_deleted_seedtissue <- vst_DE_ZE_Cache
vst_deleted_seedtissue <- vst_deleted_seedtissue[,!(vst_deleted_seedtissue != dds_deleted_seedtissue$NGI.ID)]




ndeg_Test <- sapply(9:14,function(i){
    vec <- rep(0,length(resultsNames(dds_deleted_seedtissue)))
    vec[i] <- 1
    vec[2] <- 1
    # vec[i-9] <- -1 IS REMOVING THE EFFECT OF FMG, THEREFORE WE SEE THE EXPRESSION OF ZE TISSUE AT THAT TIME
    #WE AREN'T JUST REMOVING THE TIME EFFECT, WE ARE REMOVING THE TISSUE EFFECT RELATED WITH THAT TIME INFORMATION
    vec[i-6] <- -1
#    if(i>13){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        #        vec[i-1] <- -1
        #        vec[i-10] <- -1
#    }
    print(vec)
    print(resultsNames(dds_deleted_seedtissue))
    sapply(extract_results(dds_deleted_seedtissue,vst_deleted_seedtissue,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})
















ndeg3_export_names <- resultsNames(dds_DE_ZE_Cache)[2:10]
ndeg3_export_names <- str_replace(ndeg3_export_names, "Time", "FMG_Time")
ndeg3_export_names <- str_replace(ndeg3_export_names, "vs_B1", "vs_B(-1)")

####Time vs B1 FMG Tissue
###differential expression based on Time only
ndeg3 <- sapply(3:11,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg3_export_names[i-2]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[3:11]
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
#ndeg3trans <- t(ndeg3)
#rownames(ndeg3trans) <- barnames
#colnames(ndeg3trans) <- c("all","up","dn")
#barplot(ndeg3trans,beside = TRUE, las=2, horiz=F)


#' ### Venn Diagram
#grid.newpage()
#grid.draw(venn.diagram(c("B6 vs B1"=B6vsB1[1],
#                            "B7 vs B1"=B7vsB1[1]),
#                       NULL,
#                       fill=pal[1:2]))




ndeg4_export_names <- resultsNames(dds_DE_ZE_Cache)[12:20]
ndeg4_export_names <- str_replace(ndeg4_export_names, "Time", "ZE_Time")
ndeg4_export_names <- str_replace(ndeg4_export_names, "vs_B1", "vs_B(-1)")

resultsNames(dds_DE_ZE_Cache)
####Time vs ZETissue (not adjusted for FMG vs ZE)
ndeg4 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    if(i<14){
        vec[i-9]<-1
        if(i>12){
            vec[i-10]<- -1
        }
    }else{
        vec[i]<-1
    }
#    vec[i]<-1
    if(i>13){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    print(vec)
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg4_export_names[i-11]),length)
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
#ndeg4trans <- t(ndeg4)
#rownames(ndeg4trans) <- barnames
#colnames(ndeg4trans) <- c("all","up","dn")
#barplot(ndeg4trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG



resultsNames(dds_DE_ZE_Cache)
results(dds_DE_ZE_Cache, name="TimeB2.TissueZE")
results(dds_DE_ZE_Cache, contrast=list(c("Time_B2_vs_B1","TimeB2.TissueZE")))
attr(dds_DE_ZE_Cache, "modelMatrixType")

####This should combine the main effect (Tissue vs Tissue), with the interaction effect (ZETissue.Time)
# The main effect is Tissue_ZE_vs_FMG, plus the interaction effects in the ZE tissue
# Would I still do -1 for the previous interaction term?
ndeg4devel <- sapply(12:20,function(i){
    
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[2] <- 1 #Time in 
    vec[i] <- 1
    if(i>12){
        vec[i-1] <- -1
    }
    print(vec)
    
    
    
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg4_export_names[i-11]),length)
})

barnames <- resultsNames(dds_DE_ZE_Cache)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")
barnames <- str_replace(barnames,"ZE","vs B")
barnames <- str_c(barnames[1:9],1:9)

#Differential expression ZE vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg4devel) <- barnames
rownames(ndeg4devel) <- c("all","up","dn")
barplot(ndeg4devel,beside = TRUE, las=2, horiz=F)

#To determine the Tissue Effect, without looking at Time
# Contrast Tissue, ZE and FMG
ndeg4devel_2 <- sapply(1,function(i){
    
    print(i)
    vec <- c("Tissue","ZE","FMG")
    print(vec)
    
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg4_export_names[i-11]),length)
})

#Condition effect in Tissue = Condition + Interaction term (eg TimeB2 vs TimeB1 + TissueZE.TimeB2)

















resultsNames(dds_DE_ZE_Cache_Releveled)
{
    ndeg_relevel_ZE_vs_FMG_devel <- sapply(12:20,function(i){
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
        vec[i]<-1
        #vec[i-9]<-1
        if(i>12){
            #each unit is being compared to the previous unit, with the -1 (need to adjust names)
            
            vec[i-10] <- -1
        #    vec[i-1] <- -1
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
    colnames(ndeg_relevel_ZE_vs_FMG_devel) <- barnames
    rownames(ndeg_relevel_ZE_vs_FMG_devel) <- c("all","up","dn")
    barplot(ndeg_relevel_ZE_vs_FMG_devel,beside = TRUE, las=2, horiz=F)
}

ndeg3-ndeg4
ndeg3
ndeg4


















resultsNames(dds_DE_ZE_Cache)
####Time vs ZETIssue (maybe adjusted for FMG vs ZE?)
ndeg5 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[i-9]<- -1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1
#        vec[i-10] <- -1        
    }
    print(vec)
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg5) <- resultsNames(dds_DE_ZE_Cache)[12:20]
rownames(ndeg5) <- c("all","up","dn")
barplot(ndeg5,beside = TRUE, las=2, horiz=F)


##0.5 for ZEvsFMG Tissue
ndeg6 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[11]<-0.5
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst_DE_ZE_Cache,vec,verbose = FALSE,
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
    vec[i-9]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)

        vec[i-10] <- -1
    }
    sapply(extract_results(dds_DE_ZE_Cache,vst,vec,verbose = TRUE,
                           export = FALSE,plot = FALSE, default_prefix = ndeg7_export_names[i-11]),length)
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
#ndeg7trans <- t(ndeg7)
#rownames(ndeg7trans) <- barnames
#colnames(ndeg7trans) <- c("all","up","dn")
#barplot(ndeg7trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG







##1 for ZE and FMG Time points, -1 to previous FMG and ZE time point
ndeg8 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    vec[i]<-1
    vec[i-9]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-10] <- -1     
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
{
    install_github("kassambara/factoextra")
    library("factoextra")
TT_mat <- str_c(dds_DE_ZE_Cache$Tissue,dds_DE_ZE_Cache$Time)
colnames(vst)


source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
TT_exp.specif <- expressionSpecificity(exp.mat = vst_DE_ZE_Cache, tissue = TT_mat, mode = c("local"), output = c("complete"))
    #exp.mat is vst_DE_ZE_Cache
#tissue will be combination of time and tissue

str_subset(colnames(TT_exp.specif), "n", negate = TRUE)
colnames(TT_exp.specif)

#gotpca, what lib is it from - factoextra
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
###error - 'x' must be numeric matrix in heatmap.2

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
}








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


