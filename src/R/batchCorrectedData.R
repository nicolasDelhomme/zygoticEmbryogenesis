#' ---
#' title: "Batch correction ZE and SE"
#' author: "Nicolas Delhomme && Michael Stewart"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
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
  library(parallel)
  library(plotly)
  library(tibble)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Data
load(here("analysis/salmon/ZE-29Seed-dds.rda"))
levels(dds.29z$Experiment) <- c("ZE","Seed")

#' # Normalise
vsd <- varianceStabilizingTransformation(dds.29z,blind=FALSE)

#' # Batch effect
#' ## Estimation
mat <- assay(vsd)

#' We select FMG and ZE from mature seeds (29Seed), that correspond to
#' Stage B8 and B9 of the ZE
sel <- dds.29z$Experiment == "Seed" | dds.29z$Time %in% c("B8","B9")

mat <- mat[,sel]

batch = dds.29z$Experiment[sel]

contrasts(batch) <- contr.sum(levels(batch))

design <- model.matrix(~batch)

fit <- lmFit(mat, design)

#' beta is the actual batch correction for each gene
beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0

#' cross-validation with the removeBatchEffect function
stopifnot(all((mat - beta %*% t(design[,-1,drop=FALSE]))[1:6,1:6] == 
                removeBatchEffect(mat,dds.29z$Experiment[sel])[1:6,1:6]))

#' ## Correction 
load(here("analysis/salmon/ZE-SE-dds.rda"))
levels(dds.sz$Experiment) <- c("ZE","SE")
vsd <-varianceStabilizingTransformation(dds.sz,blind=FALSE)

fullbatch <- vsd$Experiment
contrasts(fullbatch) <- contr.sum(levels(fullbatch))
fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]

assay(vsd) <- assay(vsd) - beta %*% t(fullbatch)

#' # Quality Assessment
#' ## PCA
pc <- prcomp(t(assay(vsd)))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=2

#' An the number of possible combinations
nlevel=nlevels(vsd$Tissue) * nlevels(vsd$Time)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(vsd)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Heatmap
#' 
#' Filter for noise - the filtering is meant to restrict the number
#' of genes to be visualised to a reasonable amount that can be 
#' viewed in a limited amount of time
#' 
conds <- factor(paste(vsd$Time,vsd$Tissue))
sels <- rangeFeatureSelect(counts=assay(vsd),
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 6

#' * Heatmap (>20k genes)
hm <- heatmap.2(t(scale(t(assay(vsd)[sels[[vst.cutoff+1]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,
                col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",labels=conds,cex=0.8)

#' # Export
dir.create(here("analysis/batchCorrection"),showWarnings=FALSE)
save(vsd,file=here("analysis/batchCorrection/vsd.rda"))

#' # Example limma design
#' look https://www.bioconductor.org/packages//2.11/bioc/vignettes/limma/inst/doc/usersguide.pdf
#' e.g. chapter 8.5.4
#' select the data to drop Stage B0, B9 and B10 (insufficient / lack of replication)
vsd <- vsd[,!vsd$Time %in% c("B0","B9","B10")]
vsd$Time <- droplevels(vsd$Time)
Stage<-vsd$Time
Tissue<-vsd$Tissue
design <- model.matrix(~Stage*Tissue)
fit <- lmFit(assay(vsd), design)
fit <- eBayes(fit)
colnames(design)
topTable(fit, coef=10-9)

#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
