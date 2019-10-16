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


#' ## Normalisation for visualisation
#' the normalisation is aware to take advantage of the model to determine the dispersion
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
#' * 245379 
#' The gene is medium expressed and show different patterns
#' in ECM and FLM
line_plot(dds,vst,"MA_99998g0010.1")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `ECM` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### FLM _vs._ ECM at T3


B2vsB1FMG <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                             export = FALSE,plot = FALSE,
                         default_prefix="Tissue_B2-FMG_B1-FMG__",
                         labels=paste0(colData(dds)$Tissue,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Tissue!="ZE")



####original, gives total bar plot, no distinction between up or down regulated.
ndeg <- lapply(2:3,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    if(i>2){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1 
    }
    return(list(length(extract_results(dds,vst,vec,verbose = FALSE,
                    export = FALSE,plot = FALSE)),vec))
})


barplot(unlist(ndeg),beside = TRUE,las=2)
names(ndeg) <- resultsNames(dds)[2:10]




####Time vs B1 FMG Tissue
###differential expression based on Time only
ndeg3 <- sapply(2:10,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    if(i>2){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

barnames <- resultsNames(dds)[2:10]
barnames <- str_replace_all(barnames,"_","")
barnames <- str_replace(barnames,"Time"," ")
barnames <- str_replace(barnames,"vs"," vs ")

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg3) <- barnames
rownames(ndeg3) <- c("all","up","dn")
barplot(ndeg3,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg3trans <- t(ndeg3)
rownames(ndeg3trans) <- barnames
colnames(ndeg3trans) <- c("all","up","dn")
barplot(ndeg3trans,beside = TRUE, las=2, horiz=F)

resultsNames(dds)
B6vsB1 <- sapply(extract_results(dds,vst,c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                          export = FALSE,plot = FALSE),length)
B6vsB1 <- data.frame(B6vsB1)

B7vsB1 <- sapply(extract_results(dds,vst,c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                          export = FALSE,plot = FALSE),length)
#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(c("B6 vs B1"=B6vsB1[1],
                            "B7 vs B1"=B7vsB1[1]),
                       NULL,
                       fill=pal[1:2]))






####Time vs ZETissue (not adjusted for FMG vs ZE)
ndeg4 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

barnames <- resultsNames(dds)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")

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
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-0.5
    vec[13]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg5) <- resultsNames(dds)[12:20]
rownames(ndeg5) <- c("all","up","dn")
barplot(ndeg5,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg5trans <- t(ndeg5)
rownames(ndeg5trans) <- resultsNames(dds)[12:20]
colnames(ndeg5trans) <- c("all","up","dn")
barplot(ndeg5trans,beside = TRUE, las=2, horiz=F)
    ###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG



##0.5 for ZEvsFMG Tissue
ndeg6 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    vec[11]<-0.5
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg6) <- resultsNames(dds)[12:20]
rownames(ndeg6) <- c("all","up","dn")
barplot(ndeg6,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg6trans <- t(ndeg6)
rownames(ndeg6trans) <- resultsNames(dds)[12:20]
colnames(ndeg6trans) <- c("all","up","dn")
barplot(ndeg6trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG







###comparing ZE and FMG time points, with -1 on the previous time point of FMG
ndeg7 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    vec[i-10]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)

        vec[i-11] <- -1
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

barnames <- resultsNames(dds)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")


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







##0.5 for ZE and FMG Time points, -1 to previous FMG time point
ndeg8 <- sapply(12:20,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-0.5
    vec[i-10]<-1
    if(i>12){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-11] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

barnames <- resultsNames(dds)[12:20]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,".Tissue"," ")

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











ZEB2B3 <- extract_results(dds,vst,c(0,0,0,0,0,
                                    0,0,0,0,0,
                                    1,0.5,0.5,0,0,
                                    0,0,0,0,0),
                          verbose = TRUE, export = FALSE,plot = FALSE,
                          default_prefix="Tissue_ZE_vs_S_", 
                          labels=paste0(colData(dds)$Tissue,
                                        colData(dds)$Time),
                          sample_sel=colData(dds)$Tissue!="FMG")

ZEB4B5B6 <- extract_results(dds,vst,c(0,0,0,0,0,
                                      0,0,0,0,0,
                                      1,0,0,0.5,0.5,
                                      0.5,0,0,0,0),
                            verbose = TRUE, export = FALSE,plot = FALSE,
                            default_prefix="Tissue_ZE_vs_S_", 
                            labels=paste0(colData(dds)$Tissue,
                                          colData(dds)$Time),
                            sample_sel=colData(dds)$Tissue!="FMG")

ZEB4B5B6minus <- extract_results(dds,vst,c(0,0,0,0,0,
                                      0,0,0,0,0,
                                      1,0,-1,0.5,0.5,
                                      0.5,0,0,0,0),
                            verbose = TRUE, export = FALSE,plot = FALSE,
                            default_prefix="Tissue_ZE_vs_S_", 
                            labels=paste0(colData(dds)$Tissue,
                                          colData(dds)$Time),
                            sample_sel=colData(dds)$Tissue!="FMG")


ZEB7B8B910 <- extract_results(dds,vst,c(0,0,0,0,0,
                                       0,0,0,0,0,
                                       1,0,0,0,0,
                                       0,0.5,0.5,0.5,0.5),
                              verbose = TRUE, export = FALSE,plot = FALSE,
                             default_prefix="Tissue_ZE_vs_S_", 
                             labels=paste0(colData(dds)$Tissue,
                                           colData(dds)$Time),
                             sample_sel=colData(dds)$Tissue!="FMG")

ZEB7B8B910minus <- extract_results(dds,vst,c(0,0,0,0,0,
                                        0,0,0,0,0,
                                        1,0,0,0,0,
                                        -1,0.5,0.5,0.5,0.5),
                              verbose = TRUE, export = FALSE,plot = FALSE,
                              default_prefix="Tissue_ZE_vs_S_", 
                              labels=paste0(colData(dds)$Tissue,
                                            colData(dds)$Time),
                              sample_sel=colData(dds)$Tissue!="FMG")





#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("B2-B3"=ZEB2B3$all,
                            "B4-B6"=ZEB4B5B6$all,
                            "B7-B10"=ZEB7B8B910$all),
                       NULL,
                       fill=pal[1:3]))










FMG_S <- extract_results(dds,vst,"Tissue_FMG_vs_S",
                        default_prefix="Tissue_FMG_vs_S_",
                        labels=paste0(colData(dds)$Tissue,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Tissue!="ZE")

#' ### FLM _vs._ ECM at T7
#' Here we want to conmbine the effect of FLM-ECM at time T3 and the specific
#' FLM:T7 interaction 
ZE_S <- extract_results(dds,vst,"Tissue_ZE_vs_S", # c(Tissue,"ZE","S")
                       default_prefix="Tissue_ZE_vs_S_", 
                       labels=paste0(colData(dds)$Tissue,
                                     colData(dds)$Time),
                       sample_sel=colData(dds)$Tissue!="FMG")


resultsNames(dds)






#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("1"=ZEB1B2$all,
                            "2"=ZEB4B5$all),
          NULL,
          fill=pal[1:2]))

grid.newpage()
grid.draw(venn.diagram(list("1"=ZEB4$all,
                            "2"=ZEB5$all),
                       NULL,
                       fill=pal[1:2]))


####Intercept (B1 VS B1?), B2 vs B1 and B3 VS B1
grid.newpage()
grid.draw(venn.diagram(list("B1 vs B1"=B1vsB1FMG,
                            "B2 vs B1"=B2vsB1FMG,
                            "B3 vs B1"=B3vsB1FMG),
                       NULL,
                       fill=pal[1:3]))


grid.newpage()
grid.draw(venn.diagram(list("B7 vs B1"=B7vsB1FMG,
                            "B8 vs B1"=B8vsB1FMG,
                            "B9 vs B1"=B9vsB1FMG,
                            "B10 vs B1"=B10vsB1FMG,
                            "B738 vs B1"=B738vsB1FMG),
                       NULL,
                       fill=pal[1:5]))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


