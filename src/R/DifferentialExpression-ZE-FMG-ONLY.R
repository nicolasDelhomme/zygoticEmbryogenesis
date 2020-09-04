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

##THIS IS ALL THE PREVIOUS TESTING STUFF
{
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
}
##THIS IS ALL THE PREVIOUS TESTING STUFF


#export names
ndeg_time_exnames <- resultsNames(dds_DE_ZE_Cache)[3:10]
ndeg_time_exnames <- str_replace(ndeg_time_exnames, "Time", "FMG_Time")
ndeg_time_exnames <- str_replace(ndeg_time_exnames, "vs_B1", "vs_B")
ndeg_time_exnames <- str_c(ndeg_time_exnames, 1:8)

ndeg_time_ze_exnames <- resultsNames(dds_DE_ZE_Cache)[3:10]
ndeg_time_ze_exnames <- str_replace(ndeg_time_ze_exnames, "Time", "ZE_Time")
ndeg_time_ze_exnames <- str_replace(ndeg_time_ze_exnames, "vs_B1", "vs_B")
ndeg_time_ze_exnames <- str_c(ndeg_time_ze_exnames, 1:8)

ndeg_tissue <- resultsNames(dds_DE_ZE_Cache)[2:10]
ndeg_tissue <- str_c(ndeg_tissue[1], "_B",1:9)


resultsNames(dds_DE_ZE_Cache)
resultsNames(dds_DE_ZE_Cache_Releveled)
#3 to 10 are Time vs Time
#11 to 18 are TissueZE.Time



####THIS IS RIGHT
####THIS IS RIGHT
####THIS IS RIGHT

#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_main_add_interaction <- sapply(10:18,function(i){
    vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
    if(i == 10){
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
                               export = TRUE,plot = FALSE, default_dir = here("analysis/DE/ZF_TissueEffect"), default_prefix = ndeg_tissue[i-9]),length)
})
#Barplot TISSUE ZE VS FMG - DIFFERENCE BETWEEN TISSUES
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds_DE_ZE_Cache_Releveled)[2:10]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:8],1:8)
        
        barnames <- str_c("B",1:8)
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

ndeg_time <- sapply(3:10,function(i){
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
ndeg_time_ze <- sapply(3:10,function(i){
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
{
    venn_ZF_list <- list("FMG_B4_vs_B3"=DE_FMG_B3_B4,
                         "FMG_B7_vs_B6"=DE_FMG_B7_B6,
                         "ZE_B4_vs_B3"=DE_ZE_B3_B4,
                         "ZE_B7_vs_B6"=DE_ZE_B7_B6)
    
    grid.newpage()
    png(file=here(str_c("/analysis/DE/", "venn_ZF_B3vB4_B7vB6.png")),
        width=1200, height=800)
    grid.draw(venn.diagram(venn_ZF_list,
                           NULL,
                           fill=pal[1:4]))
    dev.off()
}

venn_ZF_list <- venn(venn_ZF_list, show.plot = FALSE)
### EXTRACT THE GENES THAT ARE IN THESE VENN INTERSECTS

venn_ZF_list
colnames(venn_ZF_list)
rownames(venn_ZF_list)



venn_ZF_genes <- attr(venn_ZF_list, "intersections")
names(venn_ZF_genes)

#Enrichment of venn_ZF_genes
enr_venn_ZF_genes <- lapply(venn_ZF_genes[1:length(venn_ZF_genes)], function(x){
    print(length(x))
    if(length(x) > 1){
        x <- str_replace(x,"[.]1","")
        print(x)
        gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
    }else{
        NULL
    }
})
###SAVE VENN ZF ENR GENES
names(enr_venn_ZF_genes) <- str_replace_all(names(enr_venn_ZF_genes), ":", "_intersecting_")
names(enr_venn_ZF_genes) <- str_replace_all(names(enr_venn_ZF_genes), "B4_vs_B3", "EarlyStage")
names(enr_venn_ZF_genes) <- str_replace_all(names(enr_venn_ZF_genes), "B7_vs_B6", "LateStage")


enr_venn_ZF_genes$FMG_EarlyStage$go[1,]

enr2tsv(enr_venn_ZF_genes, file=paste0(here("/analysis/DE/ZF_Venn/"),"enrichedGenes"))

#Plot all Venn Diagram Treemaps
for(i in 1:length(enr_venn_ZF_genes)){
    
    x <- enr_venn_ZF_genes
    a <- length(x[[i]])
    dir <- "analysis/DE/ZF_Venn/"
    
    plotname <- names(x[i])
    
    if(a != 0){
        print(plotname)
        
        #plot and save go treemap
        if(is.null(nrow(x[[i]][[1]])) == FALSE){
            png(file=here(str_c(dir, "go_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
            dev.off()
            print("go")
        }
        
        #plot and save mapman treemap
        if(is.null(nrow(x[[i]][[2]])) == FALSE){
            png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
            dev.off()
            print("map")
        }
        
    }
}









ndeg_time[1,3]
ndeg_time[1,6]
ndeg_time_ze[1,3]
ndeg_time_ze[1,6]




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
TT_exp.specif <- expressionSpecificity(exp.mat = vst_DE_ZE_Cache, tissues = TT_mat, mode = "local", output = "complete")
    #exp.mat is vst_DE_ZE_Cache
#tissue will be combination of time and tissue

TT_mat_levels <- str_sort(TT_mat, numeric = TRUE)
TT_mat_levels <- unique(TT_mat_levels)

file.path(here("analysis/tissuespecificity/ZF",paste0("tissueSpecificity_ZF_",a,"_genes.csv")))
TT_exp.specif_peaks <- sapply(1:length(TT_mat_levels),function(i){
    if(i <= 3){
        a <- str_c(TT_mat_levels[i],",",TT_mat_levels[i+10])
    }else{
        a <- TT_mat_levels[i]
    }
    print(a)
    
    peakcol <- length(colnames(TT_exp.specif[,]))
    
    b <- rownames(TT_exp.specif[TT_exp.specif[,peakcol] == a,])
    if(length(b) > 0){
        write.csv(b,file=file.path(here("analysis/tissuespecificity/ZF",paste0("tissueSpecificity_ZF_",a,"_genes.csv"))))
    }
    return(b)
})

names(TT_exp.specif_peaks) <- TT_mat_levels
names(TT_exp.specif_peaks)


source(here("UPSCb-common/src/R/gopher.R"))
enr_TT_exp.specif_peaks_ZF <- lapply(TT_exp.specif_peaks[1:length(TT_exp.specif_peaks)], function(x){
    print(length(x))
    if(length(x) > 1){
        x <- str_replace(x,"[.]1","")
        print(x)
        gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
    }else{
        NULL
    }
})
###SAVE PEAK ZF ENR GENES
enr2tsv(enr_TT_exp.specif_peaks_ZF, file=paste0(here("/analysis/tissuespecificity/ZF/"),"enrichedGenes"))


#PLOTTING TREEMAPS OF SZ Tissue Specificity
for(i in 1:length(enr_TT_exp.specif_peaks_ZF)){
    
    x <- enr_TT_exp.specif_peaks_ZF
    a <- length(x[[i]])
    dir <- "analysis/tissuespecificity/ZF/Treemaps/"
    
    plotname <- names(x[i])
    
    if(a != 0){
        print(plotname)
        
        #plot and save go treemap
        if(is.null(nrow(x[[i]][[1]])) == FALSE){
            png(file=here(str_c(dir, "go_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
            dev.off()
            print("go")
        }
        
        #plot and save mapman treemap
        if(is.null(nrow(x[[i]][[2]])) == FALSE){
            png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
            dev.off()
            print("map")
        }
        
    }
}




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



#' Extraction of Genes (adapted from Elena's script)
dds_DE_ZE_Cache$Time
resultsNames(dds_DE_ZE_Cache)
{
    #' Loop which produces a variable for each time comparison, between B4 to B8, whcih contains the 'results' function result from DESeq2
    for(i in 3:11){
        
        let <- i-2
        fir <- str_c("B",(let+1))
        sec<- str_c("B",let)
        
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        vec[i] <- 1
        if(i>3){
            vec[i-1] <- -1
        }
        nam_comb <- str_c(fir, "_", sec)
        print(nam_comb)
        
        print(vec)
        
        assign(paste("res_", nam_comb, sep=""),
               as.data.frame(results(
                   dds_DE_ZE_Cache, contrast = vec,
                   filter = rowMedians(counts(dds_DE_ZE_Cache)),
                   parallel = TRUE)))
    }
    
    #' Extracts log2FoldChange and places them into a new dataframe
    zf_DESeq <- data.frame(Genes = rownames(res_B2_B1), 
                           L2FC_B2_vs_B1 = res_B2_B1$log2FoldChange,
                           L2FC_B3_vs_B2 = res_B3_B2$log2FoldChange,
                           L2FC_B4_vs_B3 = res_B4_B3$log2FoldChange,
                           L2FC_B5_vs_B4 = res_B5_B4$log2FoldChange, 
                           L2FC_B6_vs_B5 = res_B6_B5$log2FoldChange,
                           L2FC_B7_vs_B6 = res_B7_B6$log2FoldChange,
                           L2FC_B8_vs_B7 = res_B8_B7$log2FoldChange,
                           L2FC_B9_vs_B8 = res_B9_B8$log2FoldChange,
                           L2FC_B10_vs_B9 = res_B10_B9$log2FoldChange)
    zf_DESeq[is.na(zf_DESeq)] <- 0
    zf_DESeq$Genes <- sub("\\.1$", "", zf_DESeq$Genes)
    write_tsv(zf_DESeq, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_DESeq_cyto.tsv")
    
    #' Extracts padj and places them into a new dataframe
    zf_padj_df <- res_B5_B4
    zf_padj <- data.frame(Genes = rownames(res_B2_B1), 
                          padj_B2_vs_B1 = res_B2_B1$padj,
                          padj_B3_vs_B2 = res_B3_B2$padj,
                          padj_B4_vs_B3 = res_B4_B3$padj,
                          padj_B5_vs_B4 = res_B5_B4$padj,
                          padj_B6_vs_B5 = res_B6_B5$padj,
                          padj_B7_vs_B6 = res_B7_B6$padj,
                          padj_B8_vs_B7 = res_B8_B7$padj,
                          padj_B9_vs_B8 = res_B9_B8$padj,
                          padj_B10_vs_B9 = res_B10_B9$padj)
    zf_padj[is.na(zf_padj)] <- 1
    zf_padj$Genes <- sub("\\.1$", "", zf_padj$Genes)
    write.csv(zf_padj, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_padj_cyto.csv")
}
#' Use the above exported files in infomaptools

#Same as above, but for RELEVELED dds, which should give ZE vs Time
resultsNames(dds_DE_ZE_Cache_Releveled)
{
    #' Loop which produces a variable for each time comparison, between B4 to B8, whcih contains the 'results' function result from DESeq2
    for(i in 3:11){
        
        let <- i-2
        fir <- str_c("B",(let+1))
        sec<- str_c("B",let)
        
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache_Releveled)))
        vec[i] <- 1
        if(i>3){
            vec[i-1] <- -1
        }
        nam_comb <- str_c(fir, "_", sec)
        print(nam_comb)
        
        print(vec)
        
        assign(paste("res_", nam_comb, sep=""),
               as.data.frame(results(
                   dds_DE_ZE_Cache_Releveled, contrast = vec,
                   filter = rowMedians(counts(dds_DE_ZE_Cache_Releveled)),
                   parallel = TRUE)))
    }
    
    #' Extracts log2FoldChange and places them into a new dataframe
    zf_DESeq_releveled <- data.frame(Genes = rownames(res_B2_B1), 
                           L2FC_B2_vs_B1 = res_B2_B1$log2FoldChange,
                           L2FC_B3_vs_B2 = res_B3_B2$log2FoldChange,
                           L2FC_B4_vs_B3 = res_B4_B3$log2FoldChange,
                           L2FC_B5_vs_B4 = res_B5_B4$log2FoldChange, 
                           L2FC_B6_vs_B5 = res_B6_B5$log2FoldChange,
                           L2FC_B7_vs_B6 = res_B7_B6$log2FoldChange,
                           L2FC_B8_vs_B7 = res_B8_B7$log2FoldChange,
                           L2FC_B9_vs_B8 = res_B9_B8$log2FoldChange,
                           L2FC_B10_vs_B9 = res_B10_B9$log2FoldChange)
    zf_DESeq_releveled[is.na(zf_DESeq_releveled)] <- 0
    zf_DESeq_releveled$Genes <- sub("\\.1$", "", zf_DESeq_releveled$Genes)
    write_tsv(zf_DESeq_releveled, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_DESeq_releveled_cyto.tsv")
    
    #' Extracts padj and places them into a new dataframe
    zf_padj_releveled_df <- res_B5_B4
    zf_padj_releveled <- data.frame(Genes = rownames(res_B2_B1), 
                          padj_B2_vs_B1 = res_B2_B1$padj,
                          padj_B3_vs_B2 = res_B3_B2$padj,
                          padj_B4_vs_B3 = res_B4_B3$padj,
                          padj_B5_vs_B4 = res_B5_B4$padj,
                          padj_B6_vs_B5 = res_B6_B5$padj,
                          padj_B7_vs_B6 = res_B7_B6$padj,
                          padj_B8_vs_B7 = res_B8_B7$padj,
                          padj_B9_vs_B8 = res_B9_B8$padj,
                          padj_B10_vs_B9 = res_B10_B9$padj)
    zf_padj_releveled[is.na(zf_padj_releveled)] <- 1
    zf_padj_releveled$Genes <- sub("\\.1$", "", zf_padj_releveled$Genes)
    write.csv(zf_padj_releveled, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_padj_releveled_cyto.csv")
}





#' Extraction of Genes (adapted from Elena's script)
dds_DE_ZE_Cache$Time
resultsNames(dds_DE_ZE_Cache)
{
    #' Loop which produces a variable for each time comparison, between B4 to B8, whcih contains the 'results' function result from DESeq2
    for(i in 11:20){
        
        let <- i-10
        fir <- str_c("B", let)
        
        vec <- rep(0,length(resultsNames(dds_DE_ZE_Cache)))
        if(i == 11){
            vec[2] <- 1
        }else{
            vec[2] <- 1
            vec[i] <- 1
            #   vec[i-9] <- -1
        }
        nam_comb <- fir
        print(nam_comb)
        
        print(vec)
        
        assign(paste("res_", nam_comb, sep=""),
               as.data.frame(results(
                   dds_DE_ZE_Cache, contrast = vec,
                   filter = rowMedians(counts(dds_DE_ZE_Cache)),
                   parallel = TRUE)))
    }
    
    #' Extracts log2FoldChange and places them into a new dataframe
    zf_DESeq_te <- data.frame(Genes = rownames(res_B1), 
                           L2FC_B1 = res_B1$log2FoldChange,
                           L2FC_B2 = res_B2$log2FoldChange,
                           L2FC_B3 = res_B3$log2FoldChange,
                           L2FC_B4 = res_B4$log2FoldChange, 
                           L2FC_B5 = res_B5$log2FoldChange,
                           L2FC_B6 = res_B6$log2FoldChange,
                           L2FC_B7 = res_B7$log2FoldChange,
                           L2FC_B8 = res_B8$log2FoldChange,
                           L2FC_B9 = res_B9$log2FoldChange,
                           L2FC_B10 = res_B10$log2FoldChange)
    zf_DESeq_te[is.na(zf_DESeq_te)] <- 0
    zf_DESeq_te$Genes <- sub("\\.1$", "", zf_DESeq_te$Genes)
    write_tsv(zf_DESeq_te, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_DESeq_te_cyto.tsv")
    
    #' Extracts padj and places them into a new dataframe
    zf_padj_te_df <- res_B5_B4
    zf_padj_te <- data.frame(Genes = rownames(res_B1), 
                          padj_B1 = res_B1$padj,
                          padj_B2 = res_B2$padj,
                          padj_B3 = res_B3$padj,
                          padj_B4 = res_B4$padj,
                          padj_B5 = res_B5$padj,
                          padj_B6 = res_B6$padj,
                          padj_B7 = res_B7$padj,
                          padj_B8 = res_B8$padj,
                          padj_B9 = res_B9$padj,
                          padj_B10 = res_B10$padj)
    zf_padj_te[is.na(zf_padj_te)] <- 1
    zf_padj_te$Genes <- sub("\\.1$", "", zf_padj_te$Genes)
    write.csv(zf_padj_te, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/zf_padj_te_cyto.csv")
}
#' Use the above exported files in infomaptools












#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("B2-B3"=ZEB2B3$all,
                            "B4-B6"=ZEB4B5B6$all,
                            "B7-B10"=ZEB7B8B910$all),
                       NULL,
                       fill=pal[1:3]))








#treemaps
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))


enr_TT_exp.specif_peaks_ZF

#enrichment of DE genes that are in FMG vs Time
{
FL_FMG_Time <- list.files(here("analysis/DE/ZF_FMG_Time/"), 
                          recursive = TRUE, 
                          pattern = "FMG_Time",
                          full.names = TRUE)

FL_FMG_Time <- str_subset(FL_FMG_Time, "genes.csv")
#FL_FMG_Time <- str_subset(FL_FMG_Time, "up", negate = TRUE)
#FL_FMG_Time <- str_subset(FL_FMG_Time, "down", negate = TRUE)
FL_FMG_Time <- str_sort(FL_FMG_Time, numeric = TRUE)

FMG_GL <- sapply(1:length(FL_FMG_Time), function(i){
    
    if(length(read.csv(FL_FMG_Time[i])$x) == 0){
        a <- read.csv(FL_FMG_Time[i])$X
    }else{
        a <- read.csv(FL_FMG_Time[i])$x
    }
    
    a <- as.character(a)
    print(length(a))
    a <- str_replace(a,"[.]1","")
    return(a)
})
FMG_GL_Names <- sapply(1:length(FL_FMG_Time), function(i){
    a <- str_split(FL_FMG_Time, "FMG_")[[i]][3]
    a <- str_replace(a, ".csv","")
    return(a)
})
#FMG_GL_Names <- str_sort(FMG_GL_Names, numeric = TRUE)
names(FMG_GL) <- FMG_GL_Names

enr_de_fmg_time <- lapply(FMG_GL[1:length(FMG_GL)], function(x){
    print(length(x))
    if(length(x) > 1)
        gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
    else
        NULL
})
enr2tsv(enr_de_fmg_time, file=paste0(here("/analysis/DE/ZF_FMG_Time/"),"enrichedGenes"))
}



#enrichment of DE genes that are in ZE vs Time
{
    FL_ZE_Time <- list.files(here("analysis/DE/ZF_ZE_Time/"), 
                              recursive = TRUE, 
                              pattern = "ZE_Time",
                              full.names = TRUE)
    
    FL_ZE_Time <- str_subset(FL_ZE_Time, "genes.csv")
#    FL_ZE_Time <- str_subset(FL_ZE_Time, "up", negate = TRUE)
#    FL_ZE_Time <- str_subset(FL_ZE_Time, "down", negate = TRUE)
    FL_ZE_Time <- str_sort(FL_ZE_Time, numeric = TRUE)
    
    ZE_GL <- sapply(1:length(FL_ZE_Time), function(i){

        if(length(read.csv(FL_ZE_Time[i])$x) == 0){
            a <- read.csv(FL_ZE_Time[i])$X
        }else{
            a <- read.csv(FL_ZE_Time[i])$x
        }
        
        a <- as.character(a)
        a <- str_replace(a,"[.]1","")
        return(a)
    })
    ZE_GL_Names <- sapply(1:length(FL_ZE_Time), function(i){
        a <- str_split(FL_ZE_Time, "ZE_")[[i]][3]
        a <- str_replace(a, ".csv","")
        return(a)
    })
    #ZE_GL_Names <- str_sort(ZE_GL_Names, numeric = TRUE)
    names(ZE_GL) <- ZE_GL_Names
    
    enr_de_ze_time <- lapply(ZE_GL[1:length(ZE_GL)], function(x){
        print(length(x))
        if(length(x) > 1)
            gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
        else
            NULL
    })
    enr2tsv(enr_de_ze_time, file=paste0(here("/analysis/DE/ZF_ZE_Time/"),"enrichedGenes"))
}


#enrichment of DE genes that are in ZE vs FMG
{
    FL_ZE_FMG <- list.files(here("analysis/DE/ZF_TissueEffect/"), 
                             recursive = TRUE, 
                             pattern = "Tissue_ZE_vs_FMG",
                             full.names = TRUE)
    
    FL_ZE_FMG <- str_subset(FL_ZE_FMG, "genes.csv")
#    FL_ZE_FMG <- str_subset(FL_ZE_FMG, "up", negate = TRUE)
#    FL_ZE_FMG <- str_subset(FL_ZE_FMG, "down", negate = TRUE)
    FL_ZE_FMG <- str_sort(FL_ZE_FMG, numeric = TRUE)
    
    ZE_GL <- sapply(1:length(FL_ZE_FMG), function(i){
        
        if(length(read.csv(FL_ZE_FMG[i])$x) == 0){
            a <- read.csv(FL_ZE_FMG[i])$X
        }else{
            a <- read.csv(FL_ZE_FMG[i])$x
        }
        
        a <- as.character(a)
        a <- str_replace(a,"[.]1","")
        return(a)
    })
    ZE_GL_Names <- sapply(1:length(FL_ZE_FMG), function(i){
        a <- str_split(FL_ZE_FMG, "Tissue_")[[i]][2]
        a <- str_replace(a, ".csv","")
        return(a)
    })
    #ZE_GL_Names <- str_sort(ZE_GL_Names, numeric = TRUE)
    names(ZE_GL) <- ZE_GL_Names
    
    enr_de_ze_fmg <- lapply(ZE_GL[1:length(ZE_GL)], function(x){
        print(length(x))
        if(length(x) > 1)
            gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
        else
            NULL
    })
    enr2tsv(enr_de_ze_fmg, file=paste0(here("/analysis/DE/ZF_TissueEffect/"),"enrichedGenes"))
}


names(enr_de_ze_fmg)

#PLOTTING TREEMAPS OF ZF_TissueEffect
for(i in 1:length(enr_de_ze_fmg)){
    
    x <- enr_de_ze_fmg
    a <- length(x[[i]])
    dir <- "analysis/DE/ZF_TissueEffect/Treemaps/"
    
    plotname <- names(x[i])
    
    if(a != 0){
        print(plotname)
        
        #plot and save go treemap
        if(is.null(nrow(x[[i]][[1]])) == FALSE){
            png(file=here(str_c(dir, "go_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
            dev.off()
            print("go")
        }
        
        #plot and save mapman treemap
        if(is.null(nrow(x[[i]][[2]])) == FALSE){
            png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
            dev.off()
            print("map")
        }

    }
}

#PLOTTING TREEMAPS OF ZF_FMG_Time
for(i in 1:length(enr_de_fmg_time)){
    
    x <- enr_de_fmg_time
    a <- length(x[[i]])
    dir <- "analysis/DE/ZF_FMG_Time/Treemaps/"
    
    plotname <- names(x[i])
    
    if(a != 0){
        print(plotname)
        
        #plot and save go treemap
        if(is.null(nrow(x[[i]][[1]])) == FALSE){
            png(file=here(str_c(dir, "go_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
            dev.off()
            print("go")
        }
        
        #plot and save mapman treemap
        if(is.null(nrow(x[[i]][[2]])) == FALSE){
            png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
           dev.off()
            print("map")
            
        }
        
    }
}


#PLOTTING TREEMAPS OF ZF_FMG_Time
for(i in 1:length(enr_de_ze_time)){
    
    x <- enr_de_ze_time
    a <- length(x[[i]])
    dir <- "analysis/DE/ZF_ZE_Time/Treemaps/"
    
    plotname <- names(x[i])
    
    if(a != 0){
        print(plotname)
        
        #plot and save go treemap
        if(is.null(nrow(x[[i]][[1]])) == FALSE){
            png(file=here(str_c(dir, "go_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
            dev.off()
            print("go")
        }
        
        #plot and save mapman treemap
        if(is.null(nrow(x[[i]][[2]])) == FALSE){
            png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
                width=1000, height=700)
            plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
            dev.off()
            print("map")
            
        }
        
    }
}




plotEnrichedTreemap(enr_de_ze_fmg[[10]], enrichment = "go", namespace = "none")

plotname <- "testplot"




png(file=here(str_c("analysis/DE/ZF_FMG_Time/Treemaps/",plotname, ".png")),
    width=1000, height=700)
plotEnrichedTreemap(enr_de_ze_fmg[[10]], enrichment = "go", namespace = "none")
dev.off()












#plot treemaps of go and mapman for every sample

FL_ENR_FMG_Time <- list.files(here("analysis/DE/ZF_FMG_Time/"), 
                          recursive = TRUE, 
                          pattern = "enriched",
                          full.names = TRUE)

FL_ENR_FMG_Time_go <- str_subset(FL_ENR_FMG_Time, "go.tsv")
FL_ENR_FMG_Time_mapman <- str_subset(FL_ENR_FMG_Time, "mapman.tsv")

FL_ENR_FMG_Time_go <- str_sort(FL_ENR_FMG_Time_go, numeric = TRUE)
FL_ENR_FMG_Time_mapman <- str_sort(FL_ENR_FMG_Time_mapman, numeric = TRUE)

ENR_FMG_Time_Plot <- sapply(1:length(FL_ENR_FMG_Time), function(i){
    
    a <- read_tsv(FL_ENR_FMG_Time[i])
    return(a)
})













#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


