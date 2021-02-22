#' ---
#' title: "ZE and SE differential expression analysis"
#' author: "Nicolas Delhomme && Michael Stewart"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(gplots)
  library(ggplot2)
  library(here)
  library(hyperSpec)
  library(limma)
  library(readr)
  library(tibble)
  library(RColorBrewer)
  library(VennDiagram)
})

#' * Helpers
source(here("UPSCb-common/src/R/gopher.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(8,"Dark2")

#' * Data
suppressPackageStartupMessages(load(here("analysis/batchCorrection/vsd.rda")))

#' * Time points of interest
sample.sel <- (vsd$Time %in% paste0("B",3:6) & vsd$Tissue =="SE") | 
(vsd$Tissue!="SE" & vsd$Time %in% paste0("B",4:9))
vsd <- vsd[,sample.sel]
vsd$Time <- droplevels(vsd$Time)

#' * Devise a pseudotime from the PCA after bacth correction
# The order is B4-5 ZE; B3 SE; B6 ZE, B4-5 SE, B7-9 ZE, B6 SE
# corresponding to pseudo-time 1-2,3,4,5-6,7-9,10
vsd$PseudoTime <- as.integer(vsd$Time) -1 + 
  ifelse(vsd$Tissue=="SE",
         ifelse(vsd$Time=="B3",3,
                ifelse(vsd$Time %in% c("B4","B5"),4,7)),
         ifelse(vsd$Time=="B6",1,
           ifelse(vsd$Time %in% paste0("B",7:9),3,0)))
#' ```{r comment1, eval=FALSE,echo=FALSE}
#' The pseudo time calculation could be done much more accurately
#' We could check the single cell literature
#'```

#' * Gene of interests
goi <- read_tsv(here("doc/MA_list_known_SE_genes.txt"),
                col_types=cols(.default=col_character()))

#' # Functions
"line_plot" <- function(vsd=vsd,gene_id=gene_id,gene_name=gene_name){
  message(paste("Plotting",gene_id))
  sel <- grepl(gene_id,rownames(vsd))
  stopifnot(sum(sel)==1)
  
  p <- ggplot(bind_cols(as.data.frame(colData(vsd)),
                        data.frame(value=assay(vsd)[sel,])),
              aes(x=PseudoTime,y=value,col=Tissue,group=Tissue)) +
    geom_point() + geom_smooth() +
    scale_y_continuous(name="VST expression") + 
    ggtitle(label=paste("Expression for: ",gene_id,"(",gene_name,")"))
  
  suppressMessages(suppressWarnings(plot(p)))
  return(NULL)
}

#' # Plots
dev.null <- apply(goi,1,function(ro){line_plot(vsd,ro[2],ro[1])})

#' # Differential expression
#' 
#' ```{r comment2, eval=FALSE,echo=FALSE}
#' # limma 
#' # check https://www.bioconductor.org/packages//2.11/bioc/vignettes/limma/inst/doc/usersguide.pdf
#' # e.g. chapter 8.5.4
#' ```
#' We assume the model to be expression as a function of Stage and Tissue and their interaction
#' _i.e._ the two variables are not considered independent. However, we will group the Stages, to compare
#' Maturation (SE B3-5 and ZE B3-6) and Desiccation (SE B6 and ZE B7-9). 
 
vsd$Stage<-ifelse(vsd$PseudoTime<7,"Maturation","Dessication")
Stage <- vsd$Stage
Tissue<-vsd$Tissue
design <- model.matrix(~Stage*Tissue)

#' remove genes with t0o little variance
mat <- assay(vsd)[rowMads(assay(vsd))>0,]

fit <- lmFit(mat, design)
fit <- eBayes(fit)

#' The possible coefficients
colnames(design)

#' ## ZE vs. SE (maturation)
fit2 <- contrasts.fit(fit, c(0,0,0,0,-1,1))
fit2 <- eBayes(fit2)
ZSM <- topTable(fit2,p.value=0.01,lfc=0.5,adjust.method="BH",n=nrow(vsd))
sample.sel <- vsd$Stage == "Maturation" & vsd$Tissue != "FMG"
heatmap.2(t(scale(t(mat[rownames(ZSM),sample.sel]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = paste(vsd$Tissue,vsd$Time)[sample.sel],
          col=hpal,cexCol=0.8)

#' ## ZE vs. SE (desiccation)
fit2 <- contrasts.fit(fit, c(0,0,-1,1,0,0))
fit2 <- eBayes(fit2)
ZSD <- topTable(fit2,p.value=0.01,lfc=0.5,adjust.method="BH",n=nrow(vsd))
sample.sel <- vsd$Stage != "Maturation" & vsd$Tissue != "FMG"
heatmap.2(t(scale(t(mat[rownames(ZSD),sample.sel]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = paste(vsd$Tissue,vsd$Time)[sample.sel],
          col=hpal,cexCol=0.8)

#' ## Overlap
grid.newpage()
grid.draw(venn.diagram(list(ZSD=rownames(ZSD),
                            ZSM=rownames(ZSM)),NULL,fill=pal[1:2]))

#' # Gene Ontology Enrichment
GO <- lapply(lapply(list(rownames(ZSD),rownames(ZSM)),
                    sub,pattern="\\.1$",replacement=""),
             gopher,background=sub("\\.1","",rownames(mat)),
             task=list("go"),alpha=0.01,
             url="pabies")

#' # Export
dir.create(here("analysis/limma"),showWarnings=FALSE)  
write_delim(ZSD %>% rownames_to_column("ID"),here("analysis/limma/ZE-vs-SE-desiccation_DE-results.tsv"))
write_delim(ZSM %>% rownames_to_column("ID"),here("analysis/limma/ZE-vs-SE-maturation_DE-results.tsv"))
write_delim(GO[[1]]$go,here("analysis/limma/ZE-vs-SE-desiccation_DE_GO-results.tsv"))
write_delim(GO[[2]]$go,here("analysis/limma/ZE-vs-SE-maturation_DE_GO-results.tsv"))
write_delim(GO[[1]]$go[,c("id","padj")],here("analysis/limma/ZE-vs-SE-desiccation_DE_GO-for-REVIGO-results.tsv"))
write_delim(GO[[2]]$go[,c("id","padj")],here("analysis/limma/ZE-vs-SE-maturation_DE_GO-for-REVIGO-results.tsv"))

#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
