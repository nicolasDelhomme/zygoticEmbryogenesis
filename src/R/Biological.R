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

stopifnot(all(str_which(lb.filterlist, samples$NGI.ID) == 1:length(lb.filterlist)))
##assign names to the filtered filelist (which removed L002)
names(lb.filterlist) <- samples$NGI.ID
lb.filterlist <- lb.filterlist[samples$Tissue %in% c("ZE","FMG","S")] ##Are we only looking at ZE?
samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )

#' Read the expression at the gene level (there is one transcript per gene)
lb.g <- suppressMessages(tximport(files = lb.filterlist, 
                                  type = "salmon",txOut=TRUE))

counts <- round(lb.g$counts)
























#' ## Raw Data QC analysis
#' Check how many genes are never expressed - reasonable level of non-expressed genes indicated.
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
       aes(x,y,col=Tissue,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())


ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
       aes(x,y,col=Tissue,fill=)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())


#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Tissue=samples$Tissue[match(ind,samples$NGI.ID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$NGI.ID)])

#' Color by Experiment ####OBJECT TISSUE NOT FOUND############################################################
ggplot(dat,aes(x=values,group=ind,col=Tissue)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Time
ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("analysis/salmon/ZE-unnormalised-gene-expression_data.csv"))
############## change export location name

#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Time - Tissue)
##############################does not like *, + or /, will only take - as the model.
save(dds,file=here("analysis/salmon/ZE-Dataset-dds.rda"))

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is relatively variable (0 to 200 %) however the low end of the variance is likely due to poor samples where little sequencing data was actually produce that wasnt 16s bacterial.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
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

#' ### 3 first dimensions
#' Seems that different time points form small clusters and ZE and FMG tissue types appear to separate. These are Comp1 and Comp2 
#' which explains the different between most of the sampels except for one B4 ZE sample. Appears to be an outlier.
#' This seems to indicate that the Tissue and Time components explain the difference between samples.
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Time)],
              pch=c(17,19,15)[as.integer(samples$Tissue)])
legend("topleft",
       fill=pal[1:nlevels(samples$Time)],
       legend=levels(samples$Time))
legend("bottomright",
       pch=c(17,19,15),
       legend=levels(samples$Tissue))
par(mar=mar)

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

#' *  Heatmap of the 1000 most variable genes
ord <- order(rowSds(vst[sels[[vstCutoff]],]),decreasing=TRUE) [1:1000]

#' Cluster of subset compared to overall appears to be similar in shape with a shift in
#' increased VST expression and slightly more spread out density. These most variable genes
#' appear to have higher levels of VST expression.
ggplot(list(sub=rowMeans(vst[sels[[vstCutoff]],][ord,]),
            total=rowMeans(vst[sels[[vstCutoff]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' Variance in expression in different tissues can be seen. Certain tissues samples
#' appear the same as in other tissues. It appears that tissues have different
#' expression patterns. It can also be seen that over different time points, that 
#' gene expression between the samples change slightly over time.
heatmap.2(t(scale(t(vst[sels[[vstCutoff]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)


#' *  Heatmap of the 1000 least variable genes
ord <- order(rowSds(vst[sels[[vstCutoff]],])) [1:1000]

#' The subset is enriched for two subsets of genes - ones which have relatively the same
#' level of VST expression as the total gene subset but also a much higher population
#' of more highly VST expressed genes.
ggplot(list(sub=rowMeans(vst[sels[[vstCutoff]],][ord,]),
            total=rowMeans(vst[sels[[vstCutoff]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering for the least variable genes shows a very small change in separation
#' of gene expression by Tissue and Time. There is a more stark difference in the
#' Somatic samples.
heatmap.2(t(scale(t(vst[sels[[vstCutoff]],][ord,]))),
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









#import samples

#preprocess samples - separate each experiment
### Sample preprocessing
##need to delete User.ID, Sample.ID, Replicate, Mreads and X..Q30 (columns 2,4,6,7,8)
{
  samples <- filter(samples, !grepl("P11562_112",NGI.ID))
  samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )
  
  samplesZE <- subset(samples, samples$Experiment == "Zygotic Embryogenesis")
  samplesSEED <- subset(samples, samples$Experiment == "29Seed")
  samplesSE <- subset(samples, samples$Experiment == "Somatic Embryogenesis")
  samplesSEGerm <- subset(samples, samples$Experiment == "Somatic Embryogenesis Germinants")
}
#removed P11562_112 row from sample list NEED TO ADJUST HERE
samples <- filter(samples, !grepl("P11562_112",NGI.ID))
#remove this "P11562_148"



#process files using the preprocessed samples
########################ZE DATASET WITH REPLICATES ADDED AND S SPLIT INTO FMG AND ZE AT THE END
#' # Analysis
#' ## Raw data
################start of ZE
{
  lb.filelistZE <- list.files(here("data/RNA-Seq/salmon"), 
                              recursive = TRUE, 
                              pattern = "quant.sf",
                              full.names = TRUE)
  
  ##ZE is missing 204 in L002, and need to remove P11562_112 from both
  lb.filelistZE <- (str_subset(lb.filelistZE, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/RNA-Seq/salmon/P11562_112_S11_L00", negate = TRUE))
  lb.filterlistZEuniq <- (str_subset(lb.filelistZE, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/RNA-Seq/salmon/P11562_204"))
  lb.filterlistZE <- str_subset(lb.filelistZE,"L001_sortmerna_trimmomatic")
  lb.filelistZEWorking <- (str_subset(lb.filelistZE, "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/RNA-Seq/salmon/P11562_204", negate = TRUE))
  
  lb.filterlistZE1 <- str_subset(lb.filelistZEWorking,"L001_sortmerna_trimmomatic")
  names(lb.filterlistZE1) <- sapply(lapply(strsplit(lb.filterlistZE1,"_"),"[",1:2),paste,collapse="_")
  lb.filterlistZE2 <- str_subset(lb.filelistZEWorking,"L002_sortmerna_trimmomatic")
  names(lb.filterlistZE2) <- sapply(lapply(strsplit(lb.filterlistZE2,"_"),"[",1:2),paste,collapse="_")
  names(lb.filterlistZEuniq) <- sapply(lapply(strsplit(lb.filterlistZEuniq,"_"),"[",1:2),paste,collapse="_")
  
  names(lb.filterlistZE1) <- str_replace(names(lb.filterlistZE1),".*salmon/","")
  names(lb.filterlistZE2) <- str_replace(names(lb.filterlistZE2),".*salmon/","")
  names(lb.filterlistZEuniq) <- str_replace(names(lb.filterlistZEuniq),".*salmon/","")
  
  
  part1 <- suppressMessages(round(tximport(files = lb.filterlistZE1, 
                                           type = "salmon",txOut=TRUE)$counts))
  part2 <- suppressMessages(round(tximport(files = lb.filterlistZE2, 
                                           type = "salmon",txOut=TRUE)$counts))
  part3 <- suppressMessages(round(tximport(files = lb.filterlistZEuniq, 
                                           type = "salmon",txOut=TRUE)$counts)) 
  all(colnames(part1)==colnames(part2))
  countsZE <- part1 + part2
  countsREARRANGE <- cbind(countsZE[,54:55])
  countsZE <- countsZE[,-54:-55]
  countsZE <- cbind(countsZE, part3, countsREARRANGE)
  
  length(colnames(countsZE))
  length(samplesZE$NGI.ID)
  colnames(countsZE)
  samplesZE$NGI.ID
  samplesZE$NGI.ID == colnames(countsZE)
  
  samplesZE_S_to_FMG <- subset(samplesZE,samplesZE$Tissue == "S")
  samplesZE_S_to_ZE <- samplesZE_S_to_FMG
  countsZE_S_to_FMG <- countsZE
  countsZE_S_to_FMG <- countsZE_S_to_FMG[ ,samplesZE_S_to_FMG$NGI.ID]
  countsZE_S_to_ZE <- countsZE_S_to_FMG
  
  samplesZE_S_to_FMG$NGI.ID <- str_replace(samplesZE_S_to_FMG$NGI.ID,"_","_A_")
  samplesZE_S_to_FMG$Tissue <- str_replace(samplesZE_S_to_FMG$Tissue,"S","FMG")
  
  samplesZE_S_to_ZE$NGI.ID <- str_replace(samplesZE_S_to_ZE$NGI.ID,"_","_B_")
  samplesZE_S_to_ZE$Tissue <- str_replace(samplesZE_S_to_ZE$Tissue,"S","ZE")
  
  samplesZE <- subset(samplesZE, samplesZE$Tissue != "S")
  countsZE <- countsZE[ ,samplesZE$NGI.ID]
  
  colnames(countsZE_S_to_FMG) <- str_replace(colnames(countsZE_S_to_FMG),"_","_A_")
  colnames(countsZE_S_to_ZE) <- str_replace(colnames(countsZE_S_to_ZE),"_","_B_")
  
  colnames(countsZE)
  colnames(countsZE_S_to_FMG)
  colnames(countsZE_S_to_ZE)
  #append the split counts
  countsZE_Final <- cbind(countsZE,countsZE_S_to_FMG,countsZE_S_to_ZE)
  #append split samples
  samplesZE_Final <- bind_rows(samplesZE, samplesZE_S_to_FMG, samplesZE_S_to_ZE)
  
  #These final variables consist of FMG and ZE tissues.
  samplesZE_Final
  countsZE_Final
}
#################end of ZE


####start of SEED
{
  lb.filelist29Seed <- list.files("/mnt/picea/projects/spruce/uegertsdotter/29_Spruce_Seeds_Project/salmon", 
                                  recursive = TRUE, 
                                  pattern = "quant.sf",
                                  full.names = TRUE)
  
  lb.filterlist29Seed <- lb.filelist29Seed
  names(lb.filterlist29Seed) <- sapply(lapply(strsplit(lb.filterlist29Seed,"_"),"[",7:8),paste,collapse="_")
  counts29Seed <- suppressMessages(round(tximport(files = lb.filterlist29Seed, 
                                                  type = "salmon",txOut=TRUE)$counts))
  ###need to recorder the ones with a letter infront of them
  
  #removed P11562_112 row from sample list
  samples <- bind_rows(samplesZE_Final,samplesSEED)
  samples$Tissue <- factor(samples$Tissue)
  samples$Experiment <- factor(samples$Experiment)
  
  
  #####Filelists have been combined, however not able to proceed - NGI IDs do not match.
  #####All samples should have the same Tissue Type in the Somatic and 29Seed - Time is not changed either - they should all be the same.
  length(colnames(counts29Seed))
  length(colnames(countsZE_Final))
  length(samples$NGI.ID)
}
lb.g <- cbind(countsZE_Final, counts29Seed)
stopifnot(all(str_which(colnames(lb.g), samples$NGI.ID) == 1:length(colnames(lb.g))))

counts <- lb.g

################################SE IMPORTING
#####start of SE
{
  samplesSE
  
  samplesSE$Time <- samplesSE$Tissue
  samplesSE$Time <- str_c("B",samplesSE$Time)
  samplesSE$Tissue <- c("SE")
  #######PRE PROCESSING COMPLETE
  lb.filelistSomaticEmb <- list.files("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon", 
                                      recursive = TRUE, 
                                      pattern = "quant.sf",
                                      full.names = TRUE)
  lb.filterlistSomaticEmb <- str_subset(lb.filelistSomaticEmb,"L002", negate = TRUE)
  
  
  #################need to merge countsZE and part3, between column 53 and 54 (to become the new column 54)
  lb.filterlistSomaticEmb1 <- str_subset(lb.filelistSomaticEmb,"L002_sortmerna_trimmomatic")
  names(lb.filterlistSomaticEmb1) <- sapply(lapply(strsplit(lb.filterlistSomaticEmb1,"_"),"[",4:5),paste,collapse="_")
  lb.filterlistSomaticEmb2 <- str_subset(lb.filelistSomaticEmb,"L002_sortmerna_trimmomatic")
  names(lb.filterlistSomaticEmb2) <- sapply(lapply(strsplit(lb.filterlistSomaticEmb2,"_"),"[",4:5),paste,collapse="_")
  
  lb.filterlistSomaticEmb3 <- str_subset(lb.filelistSomaticEmb,"L00", negate = TRUE)
  names(lb.filterlistSomaticEmb3) <- sapply(lapply(strsplit(lb.filterlistSomaticEmb3,"_"),"[",7:8),paste,collapse="_")
  
  names(lb.filterlistSomaticEmb1) <- str_replace(names(lb.filterlistSomaticEmb1),".*salmon/","")
  names(lb.filterlistSomaticEmb2) <- str_replace(names(lb.filterlistSomaticEmb2),".*salmon/","")
  
  part1 <- suppressMessages(round(tximport(files = lb.filterlistSomaticEmb1, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  part2 <- suppressMessages(round(tximport(files = lb.filterlistSomaticEmb2, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  part3 <- suppressMessages(round(tximport(files = lb.filterlistSomaticEmb3, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  all(colnames(part1)==colnames(part2))
  countsSomaticEmb <- part1 + part2
  countsSomaticEmb <- cbind(part3, countsSomaticEmb)
}

###complete SE counts and sample processing
###merge counts and samples of ZE and SE into new set of sample and count files
{
samples_SE_plus_ZE <- bind_rows(samplesZE_Final,samplesSE)
samples_SE_plus_ZE$Tissue <- factor(samples_SE_plus_ZE$Tissue)
samples_SE_plus_ZE$Experiment <- factor(samples_SE_plus_ZE$Experiment)
samples_SE_plus_ZE$Time <- factor(samples_SE_plus_ZE$Time)

counts_SE_plus_ZE <- cbind(countsZE_Final,countsSomaticEmb)

length(colnames(counts_SE_plus_ZE))
length(samples_SE_plus_ZE$NGI.ID)

stopifnot(all(str_which(colnames(counts_SE_plus_ZE), samples_SE_plus_ZE$NGI.ID) == 1:length(colnames(counts_SE_plus_ZE))))
}

#################start of SomaticEmbGerm
{
##splitting data in order to sum the counts
lb.filterlistSomaticEmbGerm <- lb.filelistSomaticEmbGerm
lb.filterlistSomaticEmbGerm1 <- (str_subset(lb.filterlistSomaticEmbGerm, "AC7", negate = TRUE))
names(lb.filterlistSomaticEmbGerm1) <- sapply(lapply(strsplit(lb.filterlistSomaticEmbGerm1,"_"),"[",4:5),paste,collapse="_")
lb.filterlistSomaticEmbGerm2 <- (str_subset(lb.filterlistSomaticEmbGerm, "BC7", negate = TRUE))
names(lb.filterlistSomaticEmbGerm2) <- sapply(lapply(strsplit(lb.filterlistSomaticEmbGerm2,"_"),"[",4:5),paste,collapse="_")


part1 <- suppressMessages(round(tximport(files = lb.filterlistSomaticEmbGerm1, 
                                         type = "salmon",txOut=TRUE)$counts))
part2 <- suppressMessages(round(tximport(files = lb.filterlistSomaticEmbGerm2, 
                                         type = "salmon",txOut=TRUE)$counts))
all(colnames(part1)==colnames(part2))
countsEmbGerm <- part1 + part2
}
#################end of EmbGerm



#make a list of options

#counts <- ZE
#counts <- COMBINED 4
#counts <- ZE + 29SEED
#counts <- ZE + SE




#combine the necessary counts for different analysis'
###Solo ZE Dataset
###Combined 4 Datasets
###Batch correction of ZE and SE

###perform biological analysis, using standard groupings (same variable names, but rename them before every use)
{
  #' ## Raw Data QC analysis
  #' Check how many genes are never expressed - reasonable level of non-expressed genes indicated.
  sel <- rowSums(counts) == 0
  sprintf("%s%% percent (%s) of %s genes are not expressed",
          round(sum(sel) * 100/ nrow(counts),digits=1),
          sum(sel),
          nrow(counts))
  
  ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
           bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
         aes(x,y,col=Tissue,fill=Time)) + geom_col() + 
    scale_y_continuous(name="reads") +
    theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())
  
  
  ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
           bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
         aes(x,y,col=Tissue,fill=)) + geom_col() + 
    scale_y_continuous(name="reads") +
    theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())
  
  
  #' Display the per-gene mean expression
  #' 
  #' i.e. the mean raw count of every 
  #' gene across samples is calculated
  #' and displayed on a log10 scale.
  #' 
  #' The cumulative gene coverage is as expected
{
  ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
    geom_density() + ggtitle("gene mean raw counts distribution") +
    scale_x_continuous(name="mean raw counts (log10)")
  
  #' The same is done for the individual
  #' samples colored by condition. The gene coverage 
  #' across samples is extremely similar
  dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
    mutate(Tissue=samples$Tissue[match(ind,samples$NGI.ID)]) %>% 
    mutate(Time=samples$Time[match(ind,samples$NGI.ID)])
  
  #' Color by Experiment ####OBJECT TISSUE NOT FOUND############################################################
  ggplot(dat,aes(x=values,group=ind,col=Tissue)) + 
    geom_density() + ggtitle("sample raw counts distribution") +
    scale_x_continuous(name="per gene raw counts (log10)")
  
  #' Color by Time
  ggplot(dat,aes(x=values,group=ind,col=Time)) + 
    geom_density() + ggtitle("sample raw counts distribution") +
    scale_x_continuous(name="per gene raw counts (log10)")
}
}
    #' ## Export
  dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
  write.csv(counts,file=here("analysis/salmon/ZE-unnormalised-gene-expression_data.csv"))
  ############## change export location name
  
  #' ## Data normalisation 
  #' ### Preparation
  #' For visualization, the data is submitted to a variance stabilization
  #' transformation using DESeq2. The dispersion is estimated independently
  #' of the sample tissue and replicate
  #use this in (CASE)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ Time * Tissue)
  
  #use this in (CASE)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ Experiment)
  ##############################does not like *, + or /, will only take - as the model.

  #PICK ONE, DEPENDING ON WHAT IS BEING RUN
  save(dds,file=here("analysis/salmon/ZE-Dataset-dds.rda"))
  save(dds,file=here("analysis/salmon/Combined-Datasets-dds.rda"))
  save(dds,file=here("analysis/salmon/ZE-SE-Dataset-dds.rda"))
  save(dds,file=here("analysis/salmon/ZE-29Seed-dds.rda"))
  
  
  #' Check the size factors (i.e. the sequencing library size effect)
  #' 
  #' The sequencing depth is relatively variable (0 to 200 %) however the low end of the variance is likely due to poor samples where little sequencing data was actually produce that wasnt 16s bacterial.
  dds <- estimateSizeFactors(dds)
  sizes <- sizeFactors(dds)
  pander(sizes)
  boxplot(sizes, main="Sequencing libraries size factor")
  
  #' ## Variance Stabilising Transformation
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
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
  
  #' *  Heatmap of the 1000 most variable genes
  ord <- order(rowSds(vst[sels[[vstCutoff]],]),decreasing=TRUE) [1:1000]
  
  #' Cluster of subset compared to overall appears to be similar in shape with a shift in
  #' increased VST expression and slightly more spread out density. These most variable genes
  #' appear to have higher levels of VST expression.
  ggplot(list(sub=rowMeans(vst[sels[[vstCutoff]],][ord,]),
              total=rowMeans(vst[sels[[vstCutoff]],])) %>%  melt(),
         aes(x=value,col=L1)) + 
    geom_density() +
    ggtitle("Density of the subset vs. overall") + 
    scale_x_continuous(name=element_text("VST expression")) + 
    theme(legend.title=element_blank())
  
  #' Variance in expression in different tissues can be seen. Certain tissues samples
  #' appear the same as in other tissues. It appears that tissues have different
  #' expression patterns. It can also be seen that over different time points, that 
  #' gene expression between the samples change slightly over time.
  heatmap.2(t(scale(t(vst[sels[[vstCutoff]],][ord,]))),
            distfun=pearson.dist,
            hclustfun=function(X){hclust(X,method="ward.D2")},
            labRow = NA,trace = "none",
            labCol = conds,
            col=hpal)
  
  
  #' *  Heatmap of the 1000 least variable genes
  ord <- order(rowSds(vst[sels[[vstCutoff]],])) [1:1000]
  
  #' The subset is enriched for two subsets of genes - ones which have relatively the same
  #' level of VST expression as the total gene subset but also a much higher population
  #' of more highly VST expressed genes.
  ggplot(list(sub=rowMeans(vst[sels[[vstCutoff]],][ord,]),
              total=rowMeans(vst[sels[[vstCutoff]],])) %>%  melt(),
         aes(x=value,col=L1)) + 
    geom_density() +
    ggtitle("Density of the subset vs. overall") + 
    scale_x_continuous(name=element_text("VST expression")) + 
    theme(legend.title=element_blank())
  
  #' The clustering for the least variable genes shows a very small change in separation
  #' of gene expression by Tissue and Time. There is a more stark difference in the
  #' Somatic samples.
  heatmap.2(t(scale(t(vst[sels[[vstCutoff]],][ord,]))),
            distfun=pearson.dist,
            hclustfun=function(X){hclust(X,method="ward.D2")},
            labRow = NA,trace = "none",
            labCol = conds,
            col=hpal)
  

















































