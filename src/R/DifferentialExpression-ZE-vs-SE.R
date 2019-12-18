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
#' After performing batch correction (See ZE29SeedPrebatch)

#' ## Normalisation for visualisation
#' the normalisation is aware to take advantage of the model to determine the dispersion

dds.sz.ttmodel.de <- DESeqDataSetFromMatrix(
    countData = counts.sz.ttmodel,
    colData = samples.sz.ttmodel,
    design = ~Tissue * Time)
samples_SE_plus_ZE$Time
dds.sz.ttmodel.de <- dds.sz.ttmodel
dds.sz.ttmodel.de <- dds.sz.ttmodel.de[,!(dds.sz.ttmodel.de$NGI.ID == "P11562_148")]
#it appears that P464_202, P464_203, P464_205 are messed up - going to remove them. Just going to remove 205 only for now.
dds.sz.ttmodel.de <- dds.sz.ttmodel.de[,!(dds.sz.ttmodel.de$NGI.ID == "P464_205")]

colnames(dds_SE_ZE) <- dds_SE_ZE$NGI.ID


load(here("analysis/salmon/ZvsS-B4-B6-TTModel-dds.rda"))
dds.sz.ttmodel
design(dds.sz.ttmodel)

dds_SE_ZE
vst_SE_ZE

#' ## Gene of interest
#' * 245379 
#' The gene is medium expressed and show different patterns
#' in ECM and FLM
line_plot(dds_DE_ZE_Cache,vst,"MA_99998g0010.1")

#' ## Differential Expression
samples.sz.ttmodel$NGI.ID
counts.sz.ttmodel

dds.sz.ttmodel
vst.sz.ttmodel

conds <- factor(paste(dds.sz.ttmodel$Tissue,dds.sz.ttmodel$Time))
sels <- rangeFeatureSelect(counts=vst.sz.ttmodel,
                           conditions=conds,
                           nrep=3)
vstCutoff <- 5+1
vst.sz.ttmodel #66056 rows cutoff
vst.sz.ttmodel[sels[[vstCutoff]],] #14997 rows cutoff

vst.sz.ttmodel.featureselected <- vst.sz.ttmodel[sels[[vstCutoff]],]
dds.sz.ttmodel <- DESeq(dds.sz.ttmodel)


#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds.sz.ttmodel.de)

#' The model used is:
#' 
#' `Tissue * Time` meaning that the `Tissue` and `Time variable` as 
#' well as their interaction `Tissue:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model.
resultsNames(dds.sz.ttmodel.de)











#export names
ndeg_time_exnames_sz <- resultsNames(dds.sz.ttmodel)[3:6]
ndeg_time_exnames_sz <- str_replace(ndeg_time_exnames_sz, "Time", "SE_Time")
ndeg_time_exnames_sz <- str_replace(ndeg_time_exnames_sz, "vs_B4", "vs_B")
ndeg_time_exnames_sz <- str_c(ndeg_time_exnames_sz, 4:7)

ndeg_time_ze_exnames_sz <- resultsNames(dds.sz.ttmodel)[3:6]
ndeg_time_ze_exnames_sz <- str_replace(ndeg_time_ze_exnames_sz, "Time", "ZE_Time")
ndeg_time_ze_exnames_sz <- str_replace(ndeg_time_ze_exnames_sz, "vs_B4", "vs_B")
ndeg_time_ze_exnames_sz <- str_c(ndeg_time_ze_exnames_sz, 4:7)

ndeg_tissue_sz <- resultsNames(dds.sz.ttmodel)[2:6]
ndeg_tissue_sz <- str_c(ndeg_tissue_sz[1], "_B",4:8)




#ie what is the difference between genotypes ZE and FMG at the different time points
ndeg_main_SZ <- sapply(6:10,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel)))
    if(i == 6){
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
    sapply(extract_results(dds.sz.ttmodel,vst.sz.ttmodel,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_dir = here("analysis/DE/SZ_TissueEffect"), default_prefix = ndeg_tissue_sz[i-5]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds.sz.ttmodel)[7:10]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B1","vs B")
        barnames <- str_c(barnames[1:4],4:7)
        
        barnames <- str_c("B",4:8)
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_main_SZ) <- barnames
    rownames(ndeg_main_SZ) <- c("all","up","dn")
    
    barplot(ndeg_main_SZ,beside = TRUE, las=2, horiz=F, ylim = c(0,25000))
}


resultsNames(dds_DE_ZE_Cache)
resultsNames(dds.sz.ttmodel)

#############THIS IS RIGHT
#############THIS IS RIGHT
#############THIS IS RIGHT

ndeg_time_SZ <- sapply(3:6,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel)))
    vec[i] <- 1
    if(i>3){
        vec[i-1] <- -1
    }
    print(vec)
    #res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    #print(res.gene["MA_191819g0010.1",])
    
    #return(res.gene["MA_191819g0010.1",]$log2FoldChange)
    sapply(extract_results(dds.sz.ttmodel,vst.sz.ttmodel,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_dir = here("analysis/DE/SZ_SE_Time"), default_prefix = ndeg_time_exnames_sz[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds.sz.ttmodel)[3:6]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B4","vs B")
        barnames <- str_c(barnames[1:4],4:7)
    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_time_SZ) <- barnames
    rownames(ndeg_time_SZ) <- c("all","up","dn")
    
    barplot(ndeg_time_SZ,beside = TRUE, las=2, horiz=F, ylim = c(0,25000))
}

dds.sz.ttmodel_Releveled <- dds.sz.ttmodel
relevel(dds.sz.ttmodel_Releveled$Tissue, "ZE")
dds.sz.ttmodel_Releveled$Tissue <- relevel(dds.sz.ttmodel_Releveled$Tissue, "ZE")

dds.sz.ttmodel_Releveled <- DESeq(dds.sz.ttmodel_Releveled)
resultsNames(dds.sz.ttmodel_Releveled)





###Should probably relevel this one to get time and ze
ndeg_time_ze_SZ <- sapply(3:6,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel_Releveled)))
    vec[i] <- 1
    if(i>3){
        vec[i-1] <- -1
    }
    print(vec)
    #res.gene <- results(dds_DE_ZE_Cache, contrast = vec)
    #print(res.gene["MA_191819g0010.1",])
    
    #return(res.gene["MA_191819g0010.1",]$log2FoldChange)
    sapply(extract_results(dds.sz.ttmodel_Releveled,vst.sz.ttmodel,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_dir = here("analysis/DE/SZ_ZE_Time"), default_prefix = ndeg_time_ze_exnames_sz[i-2]),length)
})
#Barplot
{
    #Barplot Preprocessing
    {
        barnames <- resultsNames(dds.sz.ttmodel)[3:6]
        barnames <- str_replace_all(barnames,"_","")
        barnames <- str_replace(barnames,"Time"," ")
        barnames <- str_replace(barnames,"vs"," vs ")
        barnames <- str_replace(barnames,"vs B4","vs B")
        barnames <- str_c(barnames[1:4],4:7)

    }    
    #Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
    #differential expression based only on Time (reference is Time_B1)
    colnames(ndeg_time_ze_SZ) <- barnames
    rownames(ndeg_time_ze_SZ) <- c("all","up","dn")
    
    barplot(ndeg_time_ze_SZ,beside = TRUE, las=2, horiz=F, ylim = c(0,25000))
}

ndeg_main_SZ
ndeg_main_SZ[,3]
ndeg_time_SZ
ndeg_time_SZ[1,2]
ndeg_time_ze_SZ
ndeg_time_ze_SZ[1,2]
ndeg_time_ze_SZ_2
ndeg_time_ze_SZ_2[1,2]

###SE VS ZE TISSUE SPECIFICITY
{
    install_github("kassambara/factoextra")
    library("factoextra")
    TT_mat.sz <- str_c(dds.sz.ttmodel$Tissue,dds.sz.ttmodel$Time)
    
    
    source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
    TT_exp.specif.sz <- expressionSpecificity(exp.mat = vst.sz.ttmodel, tissues = TT_mat.sz, mode = "local", output = "complete")
    #exp.mat is vst_DE_ZE_Cache
    #tissue will be combination of time and tissue
    
    TT_mat_levels.sz <- str_sort(TT_mat.sz, numeric = TRUE)
    TT_mat_levels.sz <- unique(TT_mat_levels.sz)
    
    file.path(here("analysis/tissuespecificity/SZ",paste0("tissueSpecificity_SZ_",a,"_genes.csv")))
    TT_exp.specif_peaks.sz <- sapply(1:length(TT_mat_levels.sz),function(i){
    #    if(i <= 3){
    #        a <- str_c(TT_mat_levels[i],",",TT_mat_levels[i+10])
    #    }else{
            a <- TT_mat_levels.sz[i]
    #    }
        print(a)
        
        peakcol <- length(colnames(TT_exp.specif.sz[,]))
        
        b <- rownames(TT_exp.specif.sz[TT_exp.specif.sz[,peakcol] == a,])
        if(length(b) > 0){
            write.csv(b,file=file.path(here("analysis/tissuespecificity/SZ",paste0("tissueSpecificity_SZ_",a,"_genes.csv"))))
        }
        return(b)
    })
    
    names(TT_exp.specif_peaks.sz) <- TT_mat_levels.sz
    names(TT_exp.specif_peaks.sz)
    
    
    source(here("UPSCb-common/src/R/gopher.R"))
    enr_TT_exp.specif_peaks_SZ <- lapply(TT_exp.specif_peaks.sz[1:length(TT_exp.specif_peaks.sz)], function(x){
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
    enr2tsv(enr_TT_exp.specif_peaks_SZ, file=paste0(here("/analysis/tissuespecificity/SZ/"),"enrichedGenes"))
}

#PLOTTING TREEMAPS OF SZ Tissue Specificity
for(i in 1:length(enr_TT_exp.specif_peaks_SZ)){
    
    x <- enr_TT_exp.specif_peaks_SZ
    a <- length(x[[i]])
    dir <- "analysis/tissuespecificity/SZ/Treemaps/"
    
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














#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list("B2-B3"=ZEB2B3$all,
                            "B4-B6"=ZEB4B5B6$all,
                            "B7-B10"=ZEB7B8B910$all),
                       NULL,
                       fill=pal[1:3]))




FL_SZ <- list.files(here("analysis/DE"), 
                          recursive = TRUE, 
                          pattern = "SZ_",
                          full.names = TRUE)

FL_SZ <- str_subset(FL_SZ, "genes.csv")
FL_SZ <- str_subset(FL_SZ, "up", negate = TRUE)
FL_SZ <- str_subset(FL_SZ, "down", negate = TRUE)

FL_SZ

#SE Tissue vs TIme
FL_SZ_SE_B5_vs_B4 <- read.csv(FL_SZ[1])$X
FL_SZ_SE_B5_vs_B4 <- as.character(FL_SZ_SE_B5_vs_B4)

FL_SZ_SE_B6_vs_B5 <- read.csv(FL_SZ[2])$X
FL_SZ_SE_B6_vs_B5 <- as.character(FL_SZ_SE_B6_vs_B5)
#ZE Tissue vs Time
FL_SZ_ZE_B5_vs_B4 <- read.csv(FL_SZ[5])$X
FL_SZ_ZE_B5_vs_B4 <- as.character(FL_SZ_ZE_B5_vs_B4)

FL_SZ_ZE_B6_vs_B5 <- read.csv(FL_SZ[6])$X
FL_SZ_ZE_B6_vs_B5 <- as.character(FL_SZ_SE_B6_vs_B5)
#SE vs ZE Tissue
FL_SZ_SE.ZE_B5_vs_B4 <- read.csv(FL_SZ[3])$X
FL_SZ_SE.ZE_B5_vs_B4 <- as.character(FL_SZ_SE.ZE_B5_vs_B4)

FL_SZ_SE.ZE_B6_vs_B5 <- read.csv(FL_SZ[4])$X
FL_SZ_SE.ZE_B6_vs_B5 <- as.character(FL_SZ_SE.ZE_B6_vs_B5)

SZ_GL_Repository <- NULL
SZ_GL_Repository$Total <- list(FL_SZ_SE_B5_vs_B4,
                               FL_SZ_SE_B6_vs_B5,
                               FL_SZ_ZE_B5_vs_B4,
                               FL_SZ_SE_B6_vs_B5,
                               FL_SZ_SE.ZE_B5_vs_B4,
                               FL_SZ_SE.ZE_B6_vs_B5)

#identify unique genes expressed
#genes only in these time points within SE only
REFERENCE <- unlist(SZ_GL_Repository$Total[2:6])
SZ_SE_B5_vs_B4_GL_Uniq <- setdiff(FL_SZ_SE_B5_vs_B4,REFERENCE)
SZ_GL_Repository$SE_B5.uniq <- setdiff(FL_SZ_SE_B5_vs_B4,REFERENCE) #224 unique genes #236 uniq after only Time comparisons

REFERENCE <- unlist(c(SZ_GL_Repository$Total[1],SZ_GL_Repository$Total[3:6]))
SZ_SE_B6_vs_B5_GL_Uniq <- setdiff(FL_SZ_SE_B6_vs_B5,REFERENCE)
SZ_GL_Repository$SE_B6.uniq <- setdiff(FL_SZ_SE_B6_vs_B5,REFERENCE) #0 unique #still 0
unique(REFERENCE)
#genes only in these time points within ZE only
REFERENCE <- unlist(c(SZ_GL_Repository$Total[1:2],SZ_GL_Repository$Total[4:6]))
SZ_ZE_B5_vs_B4_GL_Uniq <- setdiff(FL_SZ_ZE_B5_vs_B4,REFERENCE)
SZ_GL_Repository$ZE_B5.uniq <- setdiff(FL_SZ_ZE_B5_vs_B4,REFERENCE) #21 unique genes #104 uniq after only Time comparisons

REFERENCE <- unlist(c(SZ_GL_Repository$Total[1:3],SZ_GL_Repository$Total[5:6]))
SZ_ZE_B6_vs_B5_GL_Uniq <- setdiff(FL_SZ_ZE_B6_vs_B5,REFERENCE)
SZ_GL_Repository$ZE_B6.uniq <- setdiff(FL_SZ_ZE_B6_vs_B5,REFERENCE) #0 unique #still 0

#genes only in these time points within SEvsZE only
REFERENCE <- unlist(c(SZ_GL_Repository$Total[1:4],SZ_GL_Repository$Total[6]))
SZ_SE.ZE_B5_vs_B4_GL_Uniq <- setdiff(FL_SZ_SE.ZE_B5_vs_B4,REFERENCE)
SZ_GL_Repository$SE.ZE_B5.uniq <- setdiff(FL_SZ_SE.ZE_B5_vs_B4,REFERENCE) #83 unique genes

REFERENCE <- unlist(c(SZ_GL_Repository$Total[1:5]))
SZ_SE.ZE_B6_vs_B5_GL_Uniq <- setdiff(FL_SZ_SE.ZE_B6_vs_B5,REFERENCE)
SZ_GL_Repository$SE.ZE_B6.uniq <- setdiff(FL_SZ_SE.ZE_B6_vs_B5,REFERENCE) #643 unique genes


#GOPHER
#PREPROCESSING
#set vst as my background gene lit
background <- rownames(vst.sz.ttmodel.featureselected)
background <- str_replace(background,"[.]1","")

#clip off .1 from the end of all genes in order to be compatible with gopher
SZ_GL_Repository$SE_B5.uniq <- str_replace(SZ_GL_Repository$SE_B5.uniq,"[.]1","")
SZ_GL_Repository$SE_B6.uniq <- str_replace(SZ_GL_Repository$SE_B6.uniq,"[.]1","")

SZ_GL_Repository$ZE_B5.uniq <- str_replace(SZ_GL_Repository$ZE_B5.uniq,"[.]1","")
SZ_GL_Repository$ZE_B6.uniq <- str_replace(SZ_GL_Repository$ZE_B6.uniq,"[.]1","")

SZ_GL_Repository$SE.ZE_B5.uniq <- str_replace(SZ_GL_Repository$SE.ZE_B5.uniq,"[.]1","")
SZ_GL_Repository$SE.ZE_B6.uniq <- str_replace(SZ_GL_Repository$SE.ZE_B6.uniq,"[.]1","")

#RUN GOPHER
#SE vs Time
GL_SE_B5B4_enr <- gopher(genes = SZ_GL_Repository$SE_B5.uniq, background = background,task = c("go","mapman"),url = "pabies")
GL_SE_B6B5_enr <- gopher(genes = SZ_GL_Repository$SE_B6.uniq, background = background,task = c("go","mapman"),url = "pabies") #invalid (0 unique)

#ZE vs Time
GL_ZE_B5B4_enr <- gopher(genes = SZ_GL_Repository$ZE_B5.uniq, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B6B5_enr <- gopher(genes = SZ_GL_Repository$ZE_B6.uniq, background = background,task = c("go","mapman"),url = "pabies") #invalid (0 unique)

#SE vs ZE
GL_SE.ZE_B5B4_enr <- gopher(genes = SZ_GL_Repository$SE.ZE_B5.uniq, background = background,task = c("go","mapman"),url = "pabies")
GL_SE.ZE_B6B5_enr <- gopher(genes = SZ_GL_Repository$SE.ZE_B6.uniq, background = background,task = c("go","mapman"),url = "pabies")

#outputs of all gophers
GL_SE_B5B4_enr
GL_SE_B6B5_enr #not valid
GL_ZE_B5B4_enr
GL_ZE_B6B5_enr #not valid
GL_SE.ZE_B5B4_enr
GL_SE.ZE_B6B5_enr


#SE_B5 GO and Mapman Treemap
plotEnrichedTreemap(GL_SE_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_SE_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B5 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#SE.ZE_B5 GO and Mapman Treemap
plotEnrichedTreemap(GL_SE.ZE_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_SE.ZE_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")
#SE.ZE_B6 GO and Mapman Treemap
plotEnrichedTreemap(GL_SE.ZE_B6B5_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_SE.ZE_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")


#' Extraction of Genes (adapted from Elena's script)
dds.sz.ttmodel$Time
{
#' Loop which produces a variable for each time comparison, between B4 to B8, whcih contains the 'results' function result from DESeq2
for(i in 4:7){
    vec1<-i+1
    vec1 <- str_c("B",vec1)
    vec2<-i
    vec2<- str_c("B",vec2)
    
    assign(paste("res_", str_c(vec1,"_",vec2), sep=""),
           as.data.frame(results(
               dds.sz.ttmodel, contrast = c("Time", vec1, vec2),
               filter = rowMedians(counts(dds.sz.ttmodel)),
               parallel = TRUE)))
}

#' Extracts log2FoldChange and places them into a new dataframe
sz_DESeq <- data.frame(Genes = rownames(res_B5_B4), 
                         L2FC_B5_vs_B4 = res_B5_B4$log2FoldChange, 
                         L2FC_B6_vs_B5 = res_B6_B5$log2FoldChange,
                         L2FC_B7_vs_B6 = res_B7_B6$log2FoldChange,
                         L2FC_B8_vs_B7 = res_B8_B7$log2FoldChange)
sz_DESeq[is.na(sz_DESeq)] <- 0
sz_DESeq$Genes <- sub("\\.1$", "", sz_DESeq$Genes)
write_tsv(sz_DESeq, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/sz_DESeq_cyto.tsv")

#' Extracts padj and places them into a new dataframe
sz_padj_df <- res_B5_B4
sz_padj <- data.frame(Genes = rownames(res_B5_B4), 
                      padj_B5_vs_B4 = res_B5_B4$padj, 
                      padj_B6_vs_B5 = res_B6_B5$padj,
                      padj_B7_vs_B6 = res_B7_B6$padj,
                      padj_B8_vs_B7 = res_B8_B7$padj)
sz_padj[is.na(sz_padj)] <- 1
sz_padj$Genes <- sub("\\.1$", "", sz_padj$Genes)
write.csv(sz_padj, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/infomap/sz_padj_cyto.csv")
}
#' Use the above exported files in infomaptools

















#treemaps
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))

#enrichment of DE genes that are in SE vs Time
{
    FL_SE_Time <- list.files(here("analysis/DE/SZ_SE_Time/"), 
                              recursive = TRUE, 
                              pattern = "SE_Time",
                              full.names = TRUE)
    
    FL_SE_Time <- str_subset(FL_SE_Time, "genes.csv")
#    FL_SE_Time <- str_subset(FL_SE_Time, "up", negate = TRUE)
#    FL_SE_Time <- str_subset(FL_SE_Time, "down", negate = TRUE)
    FL_SE_Time <- str_sort(FL_SE_Time, numeric = TRUE)
    
    SE_GL <- sapply(1:length(FL_SE_Time), function(i){
        
        if(length(read.csv(FL_SE_Time[i])$x) == 0){
            a <- read.csv(FL_SE_Time[i])$X
        }else{
            a <- read.csv(FL_SE_Time[i])$x
        }
        
        a <- as.character(a)
        a <- str_replace(a,"[.]1","")
        return(a)
    })
    SE_GL_Names <- sapply(1:length(FL_SE_Time), function(i){
        a <- str_split(FL_SE_Time, "SE_")[[i]][3]
        a <- str_replace(a, ".csv","")
        return(a)
    })
    #SE_GL_Names <- str_sort(SE_GL_Names, numeric = TRUE)
    names(SE_GL) <- SE_GL_Names
    
    enr_de_se_time <- lapply(SE_GL[1:length(SE_GL)], function(x){
        print(length(x))
        if(length(x) > 1)
            gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
        else
            NULL
    })
    enr2tsv(enr_de_se_time, file=paste0(here("/analysis/DE/SZ_SE_Time/"),"enrichedGenes"))
}

#enrichment of DE genes that are in ZE vs Time
{
    FL_ZE_Time <- list.files(here("analysis/DE/SZ_ZE_Time/"), 
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
    enr2tsv(enr_de_ze_time, file=paste0(here("/analysis/DE/SZ_ZE_Time/"),"enrichedGenes"))
}


#enrichment of DE genes that are in ZE vs SE
{
    FL_ZE_SE <- list.files(here("analysis/DE/SZ_TissueEffect/"), 
                            recursive = TRUE, 
                            pattern = "Tissue_ZE_vs_SE",
                            full.names = TRUE)
    
    FL_ZE_SE <- str_subset(FL_ZE_SE, "genes.csv")
#    FL_ZE_SE <- str_subset(FL_ZE_SE, "up", negate = TRUE)
#    FL_ZE_SE <- str_subset(FL_ZE_SE, "down", negate = TRUE)
    FL_ZE_SE <- str_sort(FL_ZE_SE, numeric = TRUE)
    
    ZE_GL <- sapply(1:length(FL_ZE_SE), function(i){
        
        if(length(read.csv(FL_ZE_SE[i])$x) == 0){
            a <- read.csv(FL_ZE_SE[i])$X
        }else{
            a <- read.csv(FL_ZE_SE[i])$x
        }
        
        a <- as.character(a)
        a <- str_replace(a,"[.]1","")
        return(a)
    })
    ZE_GL_Names <- sapply(1:length(FL_ZE_SE), function(i){
        a <- str_split(FL_ZE_SE, "Tissue_")[[i]][2]
        a <- str_replace(a, ".csv","")
        return(a)
    })
    #ZE_GL_Names <- str_sort(ZE_GL_Names, numeric = TRUE)
    names(ZE_GL) <- ZE_GL_Names
    
    enr_de_ze_se <- lapply(ZE_GL[1:length(ZE_GL)], function(x){
        print(length(x))
        if(length(x) > 1)
            gopher(x, task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
        else
            NULL
    })
    enr2tsv(enr_de_ze_se, file=paste0(here("/analysis/DE/SZ_TissueEffect/"),"enrichedGenes"))
}






#VENN DIAGRAM
{
DE_SE_B6_B5 <- list.files(here("analysis/DE/SZ_SE_Time"), 
                          recursive = TRUE, 
                          pattern = "",
                          full.names = TRUE)

DE_SE_B6_B5 <- str_subset(DE_SE_B6_B5, "Time_B6_vs_B5genes.csv")
DE_SE_B6_B5 <- sapply(1:length(DE_SE_B6_B5),function(i){
    
    GeneL <- read.csv(DE_SE_B6_B5[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})


DE_SE_B7_B6 <- list.files(here("analysis/DE/SZ_SE_Time"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)

DE_SE_B7_B6 <- str_subset(DE_SE_B7_B6, "Time_B7_vs_B6genes.csv")
DE_SE_B7_B6 <- sapply(1:length(DE_SE_B7_B6),function(i){
    
    GeneL <- read.csv(DE_SE_B7_B6[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})


DE_ZE_B6_B5 <- list.files(here("analysis/DE/SZ_ZE_Time"), 
                          recursive = TRUE, 
                          pattern = "",
                          full.names = TRUE)

DE_ZE_B6_B5 <- str_subset(DE_ZE_B6_B5, "Time_B6_vs_B5genes.csv")
DE_ZE_B6_B5 <- sapply(1:length(DE_ZE_B6_B5),function(i){
    
    GeneL <- read.csv(DE_ZE_B6_B5[i])$X
    GeneL <- as.character(GeneL)
    GeneL <- str_replace(GeneL,"[.]1","")
    
    print(GeneL)
    return(GeneL)
    
})
#ZE B6

DE_ZE_B7_B6 <- list.files(here("analysis/DE/SZ_ZE_Time"), 
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
grid.draw(venn.diagram(list("SE_B6_vs_B5"=DE_SE_B6_B5,
                            "SE_B6_vs_B7"=DE_SE_B7_B6,
                            "ZE_B6_vs_B5"=DE_ZE_B6_B5,
                            "ZE_B6_vs_B7"=DE_ZE_B7_B6),
                       
                       NULL,
                       fill=pal[1:4]))
}

#PLOTTING TREEMAPS OF SZ_TissueEffect
for(i in 1:length(enr_de_ze_se)){
    
    x <- enr_de_ze_se
    a <- length(x[[i]])
    dir <- "analysis/DE/SZ_TissueEffect/Treemaps/"
    
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

#PLOTTING TREEMAPS OF SZ_SE_Time
for(i in 1:length(enr_de_se_time)){
    
    x <- enr_de_se_time
    a <- length(x[[i]])
    dir <- "analysis/DE/SZ_SE_Time/Treemaps/"
    
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

#PLOTTING TREEMAPS OF SZ_ZE_Time
for(i in 1:length(enr_de_ze_time)){
    
    x <- enr_de_ze_time
    a <- length(x[[i]])
    dir <- "analysis/DE/SZ_ZE_Time/Treemaps/"
    
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






#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


