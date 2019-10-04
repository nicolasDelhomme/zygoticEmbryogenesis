#' ---
#' title: "T89 and _Laccaria bicolor_ Biological QA - January data"
#' author: "Nicolas Delhomme && Iryna Shutava"
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
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(here))

#' * Helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
pal <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information ########### need sample info?
samples <- read_csv("~/Git/zygoticEmbryogenesis/doc/ZE_Dataset_v3.csv",
                    col_types = cols(col_character(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_double(),
                                     col_double())) %>% 
  mutate(Tissue=factor(Tissue)) %>% 
  mutate(Time=factor(Time))

#' # Analysis
#' ## Raw data
lb.filelist <- list.files(here("data/RNA-Seq/salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Filter and select only one duplicate from each sample.
lb.filterlist <- str_subset(lb.filelist,"L001_sortmerna_trimmomatic")

#removed P11562_112 trimmomatic dataset
lb.filterlist <- (str_subset(lb.filterlist, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/RNA-Seq/salmon/P11562_112_S11_L001_sortmerna_trimmomatic/quant.sf", negate = TRUE))
#removed P11562_112 row from sample list
samples <- filter(samples, !grepl("112",NGI.ID))
#samples <- filter(samples, !grepl("P11562_201",NGI.ID)) NOT REMOVING THIS SAMPLE ANYMORE, FIXED FROM FMG TO S
#lb.filterlist <- (str_subset(lb.filterlist, "P11562_201", negate = TRUE)) NOT REMOVING THIS SAMPLE ANYMORE, FIXED FROM FMG TO S




##extract S tissue rows
##remove those rows from original
trim_samples <- subset(samples, samples$Tissue == "S")
samples <- subset(samples, samples$Tissue != "S")
samples
#trim_samples

##match NGI.IDs to filelist
##remove files
lb.filterlist_trim <- lb.filterlist
lb.filterlist_trim <- str_subset(lb.filterlist_trim, c(trim_samples$NGI.ID))
lb.filterlist_trim <- append(lb.filterlist_trim, c(str_subset(lb.filterlist, "P11562_201")))
#
#
lb.filterlist <- str_subset(lb.filterlist, c(trim_samples$NGI.ID), negate = TRUE)
lb.filterlist <- str_subset(lb.filterlist, "P11562_201", negate = TRUE)
#
#
#double_A <- trim_samples
#double_B <- trim_samples
#double_A$NGI.ID <- str_replace(double_A$NGI.ID,"P11562_","P11562_A_")
#double_B$NGI.ID <- str_replace(double_B$NGI.ID,"P11562_","P11562_B_")
#double_A$Tissue <- str_replace(double_A$Tissue,"S","FMG")
#double_B$Tissue <- str_replace(double_B$Tissue,"S","ZE")
#double_A$Replicate <- str_replace(double_A$Replicate,"S","FMG")
#double_B$Replicate <- str_replace(double_B$Replicate,"S","ZE")
#
#lb.filterlist_trim_A <- lb.filterlist_trim
#lb.filterlist_trim_B <- lb.filterlist_trim







stopifnot(all(str_which(lb.filterlist, samples$NGI.ID) == 1:length(lb.filterlist)))
##assign names to the filtered filelist (which removed L002)
names(lb.filterlist) <- samples$NGI.ID
#names(lb.filterlist_trim_A) <- double_A$NGI.ID
#names(lb.filterlist_trim_B) <- double_B$NGI.ID

lb.filterlist <- lb.filterlist[samples$Tissue %in% c("ZE","FMG")] ##Are we only looking at ZE?
#lb.filterlist_trim_A <- lb.filterlist_trim_A[double_A$Tissue %in% c("FMG")] ##Are we only looking at ZE?
#lb.filterlist_trim_B <- lb.filterlist_trim_B[double_B$Tissue %in% c("ZE")] ##Are we only looking at ZE?


#samples <- bind_rows(samples,double_A,double_B)
samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )

#, lb.filterlist_trim_A, lb.filterlist_trim_B


#' Read the expression at the gene level (there is one transcript per gene)
lb.g <- suppressMessages(tximport(files = c(lb.filterlist), 
                                  type = "salmon",txOut=TRUE))

counts <- round(lb.g$counts)

#' ## Raw Data QC analysis
#' Check how many genes are never expressed - reasonable level of non-expressed genes indicated.
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Tissue=samples$Tissue[match(ind,samples$NGI.ID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$NGI.ID)])

#' ## Export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("analysis/salmon/ZE-ZF-unnormalised-gene-expression_data.csv"))
############## change export location name

#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Time * Tissue)

#want to duplicate S here, rename each twin to FMG and ZE, and give them unique NGI.IDs



save(dds,file=here("analysis/salmon/ZE-ZF-Dataset-dds.rda"))

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is relatively variable (0 to 200 %) however the low end of the variance is likely due to poor samples where little sequencing data was actually produce that wasnt 16s bacterial.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked very well, the data is nearly homoscesdastic
#' 
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))

percent <- round(summary(pc)$importance[2,]*100)
plot(cumsum(percent))

#' ### 3 first dimensions
#' Seems that different time points form small clusters and ZE and FMG tissue types appear to separate. These are Comp1 and Comp2 
#' which explains the different between most of the sampels except for one B4 ZE sample. Appears to be an outlier.
#' This seems to indicate that the Tissue and Time components explain the difference between samples.
#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))

interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                       yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))



#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds <- factor(paste(samples$Tissue,samples$Time))
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vstCutoff <- 4+1

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples.
#' It appears that generally there is little difference between the samples across all genes
#'  - a small difference is noticable between ZE tissues and FMG-S tissues, however 
#' generally the expression levels are relatively balanced apart from the one sample 
#' in FMG B4. Somatic samples had a small section of more highly expressed genes 
#' compared to ZE and FMG.
heatmap.2(t(scale(t(vst[sels[[vstCutoff]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)


#' ## Conclusion ###################### make separate notes for this.
#' The data appears to show a correlation between Tissue and Time, with some clustering in the PCA plot
#' showing a movement of time in ascending order from right to left (+ to -) along the X axis. There can also
#' be seen that along the Y axis, there are subgroups which correlate to the Tissue type.
#' Along with this, it appears that overall gene expression does not appear starkly different but
#' slight differences can be observed between Tissue types. It appears that the most difference between
#' Tissues can be seen in the 1000 most variable genes, with some small changes between Time points within
#' those Tissue groups.
#' The final heat map of the 1000 least variable genes appear to show a very mixed expression pattern, similar
#' in appearance to the total gene heat map, between Tissue groups and Time points, except for in 
#' the Somatic tissue type.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
