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

conds <- factor(paste(dds.sz.ttmodel.de$Tissue,dds.sz.ttmodel.de$Time))
sels <- rangeFeatureSelect(counts=vst.sz.ttmodel,
                           conditions=conds,
                           nrep=3)
vstCutoff <- 6+1
vst.sz.ttmodel #66056 rows cutoff
vst.sz.ttmodel[sels[[vstCutoff]],] #11609 rows cutoff

vst.sz.ttmodel.featureselected <- vst.sz.ttmodel[sels[[vstCutoff]],]
dds.sz.ttmodel.de <- DESeq(dds.sz.ttmodel.de)


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






ndeg3.se_export_names <- resultsNames(dds.sz.ttmodel.de)[3:4]
ndeg3.se_export_names <- str_replace(ndeg3.se_export_names, "Time", "SZ_SE_Time")
ndeg3.se_export_names <- str_replace(ndeg3.se_export_names, "vs_B4", "vs_B(-1)")


#modify code - need to make sure analysis is based on the previous time (eg, look at B4-B6 range, therefore B5vsB4, B6vsB5)
#SE vs Time, is 3:4, -1 after 3
#ZE vs Time is 7:8, -1 after 7
#SE vs ZE is 7:8, 3:4, and -1 after 3


####Time vs B1 SE Tissue
###differential expression based on Time only
ndeg3.se <- sapply(3:4,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel.de)))
    vec[i]<-1
    if(i>3){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds.sz.ttmodel.de,vst.sz.ttmodel.featureselected,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg3.se_export_names[i-2]),length)
})

barnames <- resultsNames(dds.sz.ttmodel.de)[3:4]
barnames <- str_replace_all(barnames,"_","")
barnames <- str_replace(barnames,"Time"," ")
barnames <- str_replace(barnames,"vs"," vs ")
barnames <- str_replace(barnames,"vs B4","vs B")
barnames <- str_c(barnames[1:2],4:5)

#Differential expression FMG vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg3.se) <- barnames
rownames(ndeg3.se) <- c("all","up","dn")
barplot(ndeg3.se,beside = TRUE, las=2, horiz=F)


#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(c("B6 vs B1"=B6vsB1[1],
                            "B7 vs B1"=B7vsB1[1]),
                       NULL,
                       fill=pal[1:2]))




ndeg4.se_export_names <- resultsNames(dds.sz.ttmodel.de)[3:4]
ndeg4.se_export_names <- str_replace(ndeg4.se_export_names, "Time", "SZ_ZE_Time")
ndeg4.se_export_names <- str_replace(ndeg4.se_export_names, "vs_B4", "vs_B(-1)")


####Time vs ZETissue (not adjusted for FMG vs ZE)
ndeg4.se <- sapply(7:8,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel.de)))
    vec[i]<-1
    if(i>7){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)
        vec[i-1] <- -1        
    }
    sapply(extract_results(dds.sz.ttmodel.de,vst.sz.ttmodel.featureselected,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg4.se_export_names[i-6]),length)
})


barnames <- resultsNames(dds.sz.ttmodel.de)[7:8]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,"TissueZE.","")
barnames <- str_c(barnames[1:2]," vs B")
barnames <- str_c(barnames[1:2],4:5)

#Differential expression ZE vs Time (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...)
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg4.se) <- barnames
rownames(ndeg4.se) <- c("all","up","dn")
barplot(ndeg4.se,beside = TRUE, las=2, horiz=F)


ndeg7.se_export_names <- resultsNames(dds.sz.ttmodel.de)[3:4]
ndeg7.se_export_names <- str_replace(ndeg7.se_export_names, "Time", "SZ_SE_vs_ZE")
ndeg7.se_export_names <- str_replace(ndeg7.se_export_names, "vs_B4", "vs_B(-1)")

###comparing ZE and FMG time points, with -1 on the previous time point of FMG
ndeg7.se <- sapply(7:8,function(i){
    vec <- rep(0,length(resultsNames(dds.sz.ttmodel.de)))
    vec[i]<-1
    vec[i-4]<-1
    if(i>7){
        #each unit is being compared to the previous unit, with the -1 (need to adjust names)

        vec[i-5] <- -1
    }
    sapply(extract_results(dds.sz.ttmodel.de,vst.sz.ttmodel.featureselected,vec,verbose = TRUE,
                           export = TRUE,plot = FALSE, default_prefix = ndeg7.se_export_names[i-6]),length)
})

barnames <- resultsNames(dds.sz.ttmodel.de)[5:6]
barnames <- str_replace(barnames,"Time","")
barnames <- str_replace(barnames,"TissueZE.","")
barnames <- str_c(barnames[1:2]," vs B")
barnames <- str_c(barnames[1:2],4:5)

#Differential expression ZE vs FMG (redo the figure, B2vsB1, B3vsB2, B4vsB3, etc...) - looking at the tissue effect itself
#the tissue differences between FMG and ZE
#differential expression based only on Time (reference is Time_B1)
colnames(ndeg7.se) <- barnames
rownames(ndeg7.se) <- c("all","up","dn")
barplot(ndeg7.se,beside = TRUE, las=2, horiz=F)

combinedrows.se <- cbind(ndeg3.se,ndeg4.se,ndeg7.se)
barplot(combinedrows.se,beside = TRUE, las=2, horiz=F)


ndeg3.se
ndeg4.se
ndeg7.se



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
dds.sz.ttmodel.de$Time
{
#' Loop which produces a variable for each time comparison, between B4 to B8, whcih contains the 'results' function result from DESeq2
for(i in 4:7){
    vec1<-i+1
    vec1 <- str_c("B",vec1)
    vec2<-i
    vec2<- str_c("B",vec2)
    
    assign(paste("res_", str_c(vec1,"_",vec2), sep=""),
           as.data.frame(results(
               dds.sz.ttmodel.de, contrast = c("Time", vec1, vec2),
               filter = rowMedians(counts(dds.sz.ttmodel.de)),
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








#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


