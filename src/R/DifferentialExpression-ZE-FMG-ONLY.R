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
##############change file directory?

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

#' * 676331
#' The gene is very lowly expressed with some samples having no expression
line_plot(dds,vst,"676331")

#' * 315313
#' Same as 676331
line_plot(dds,vst,"Lacbi1.eu2.Lbscf0069g00950")

#' * 315258
#' The gene has a low expression and there is overrlap between the two
#' experiments, but there is visible expression difference
line_plot(dds,vst,"315258")

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


B1vsB1FMG <- extract_results(dds,vst,c(1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                             export = FALSE,plot = FALSE,
                             default_prefix="Tissue_B1-FMG_B1-FMG__",
                             labels=paste0(colData(dds)$Tissue,
                                           colData(dds)$Time),
                             sample_sel=colData(dds)$Tissue!="ZE")

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
        vec[i-1] <- -1 
    }
    return(list(length(extract_results(dds,vst,vec,verbose = FALSE,
                    export = FALSE,plot = FALSE)),vec))
})


barplot(unlist(ndeg),beside = TRUE,las=2)
names(ndeg) <- resultsNames(dds)[2:12]


###differential expression based on Time only
ndeg3 <- sapply(2:12,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    if(i>2){
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg3) <- resultsNames(dds)[2:12]
rownames(ndeg3) <- c("all","up","dn")
barplot(ndeg3,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg3trans <- t(ndeg3)
rownames(ndeg3trans) <- resultsNames(dds)[2:12]
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







ndeg4 <- sapply(14:23,function(i){
    vec <- rep(0,length(resultsNames(dds)))
    vec[i]<-1
    if(i>14){
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds,vst,vec,verbose = FALSE,
                           export = FALSE,plot = FALSE),length)
})

#differential expression based only on Time (reference is Time_B1)
colnames(ndeg4) <- resultsNames(dds)[14:23]
rownames(ndeg4) <- c("all","up","dn")
barplot(ndeg4,beside = TRUE, las=2, horiz=F)

#same as above, but showing in groups of "all" "up" and "down" regulated.
ndeg4trans <- t(ndeg4)
rownames(ndeg4trans) <- resultsNames(dds)[14:23]
colnames(ndeg4trans) <- c("all","up","dn")
barplot(ndeg4trans,beside = TRUE, las=2, horiz=F)
###appears to only have DE between B4 to B7 in ZE Tissue vs B1FMG












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


FMG_ZE <- extract_results(dds,vst,list("Tissue_ZE_vs_S","Tissue_FMG_vs_S"), #c(0,-1,1,0,0,0,0,0....)
                          default_prefix="Tissue_ZE_vs_FMG_",
                          labels=paste0(colData(dds)$Tissue,
                                        colData(dds)$Time),
                          sample_sel=colData(dds)$Tissue!="S")

FMG_ZE <- extract_results(dds,vst,c("Tissue","ZE","FMG"), #c(0,-1,1,0,0,0,0,0....)
                          export = FALSE,plot = FALSE,
                          default_prefix="Tissue_ZE_vs_FMG_",
                          labels=paste0(colData(dds)$Tissue,
                                        colData(dds)$Time),
                          sample_sel=colData(dds)$Tissue!="S")




#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("B2 vs B1"=B2vsB1FMG,
                            "B3 vs B1"=B3vsB1FMG,
                            "B4 vs B1"=B4vsB1FMG,
                            "B5 vs B1"=B5vsB1FMG,
                            "B6 vs B1"=B6vsB1FMG),
          NULL,
          fill=pal[1:5]))


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


