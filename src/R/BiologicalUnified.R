#' title: "Unified Biological QA Script"
#' author: "Michael Stwewart"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
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
samples <- read_csv("~/Git/zygoticEmbryogenesis/doc/4Datasets_v6.csv",
                    col_types = cols(col_character(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_double(),
                                     col_double(),
                                     col_factor())) %>% 
  mutate(Tissue=factor(Tissue)) %>% 
  mutate(Time=factor(Time)) %>%
  mutate(Experiment=factor(Experiment))

#' Sample preprocessing
##need to delete User.ID, Sample.ID, Replicate, Mreads and X..Q30 (columns 2,4,6,7,8)
samples <- filter(samples, !grepl("P11562_112",NGI.ID))
samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )

samplesZE <- subset(samples, samples$Experiment == "Zygotic Embryogenesis")
samples29SEED <- subset(samples, samples$Experiment == "29Seed")
samplesSE <- subset(samples, samples$Experiment == "Somatic Embryogenesis")
samplesSE_Germ <- subset(samples, samples$Experiment == "Somatic Embryogenesis Germinants")

#remove these:
#P11562_148
#P464_205



#Trim the data (removing )



#ZE Dataset (FMG and ZE Tissue)
{
  ########################ZE DATASET WITH REPLICATES ADDED AND S SPLIT INTO FMG AND ZE AT THE END
  #' # Analysis
  #' ## Raw data
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
  
  #remove "S" Tissue
  samplesZE <- subset(samplesZE, samplesZE$Tissue != "S")
  countsZE <- countsZE[ ,samplesZE$NGI.ID]
  
  colnames(countsZE_S_to_FMG) <- str_replace(colnames(countsZE_S_to_FMG),"_","_A_")
  colnames(countsZE_S_to_ZE) <- str_replace(colnames(countsZE_S_to_ZE),"_","_B_")
  
  colnames(countsZE)
  colnames(countsZE_S_to_FMG)
  colnames(countsZE_S_to_ZE)
  #append the split counts
  countsZE_Final <- cbind(countsZE,countsZE_S_to_FMG,countsZE_S_to_ZE)
  colnames(countsZE_Final)
  countsZE_Final
  #append split samples
  samplesZE_Final <- bind_rows(samplesZE, samplesZE_S_to_FMG, samplesZE_S_to_ZE)
  samplesZE_Final
  #################end of ZE
}

#29Seed Dataset
{
  lb.filelist29Seed <- list.files("/mnt/picea/projects/spruce/uegertsdotter/29_Spruce_Seeds_Project/salmon", 
                                  recursive = TRUE, 
                                  pattern = "quant.sf",
                                  full.names = TRUE)
  
  lb.filterlist29Seed <- lb.filelist29Seed
  names(lb.filterlist29Seed) <- sapply(lapply(strsplit(lb.filterlist29Seed,"_"),"[",7:8),paste,collapse="_")
  counts29Seed <- suppressMessages(round(tximport(files = lb.filterlist29Seed, 
                                                  type = "salmon",txOut=TRUE)$counts))
}

#SE Dataset
{
  #SE Dataset
  {
  ################################SE IMPORTING
  #' * Metadata
  #' Sample information ########### need sample info?
  samplesSE
  
  samplesSE$Time <- samplesSE$Tissue
  samplesSE$Time <- str_c("B",samplesSE$Time)
  samplesSE$Tissue <- c("SE")
  #######PRE PROCESSING COMPLETE
  lb.lb.filelistSE <- list.files("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon", 
                                      recursive = TRUE, 
                                      pattern = "quant.sf",
                                      full.names = TRUE)
  lb.filterlistSomaticEmb <- str_subset(lb.lb.filelistSE,"L002", negate = TRUE)
  
  
  #################need to merge countsZE and part3, between column 53 and 54 (to become the new column 54)
  lb.filterlistSE_1 <- str_subset(lb.lb.filelistSE,"L002_sortmerna_trimmomatic")
  names(lb.filterlistSE_1) <- sapply(lapply(strsplit(lb.filterlistSE_1,"_"),"[",4:5),paste,collapse="_")
  lb.filterlistSE_2 <- str_subset(lb.lb.filelistSE,"L002_sortmerna_trimmomatic")
  names(lb.filterlistSE_2) <- sapply(lapply(strsplit(lb.filterlistSE_2,"_"),"[",4:5),paste,collapse="_")
  
  lb.filterlistSE_3 <- str_subset(lb.lb.filelistSE,"L00", negate = TRUE)
  names(lb.filterlistSE_3) <- sapply(lapply(strsplit(lb.filterlistSE_3,"_"),"[",7:8),paste,collapse="_")
  
  names(lb.filterlistSE_1) <- str_replace(names(lb.filterlistSE_1),".*salmon/","")
  names(lb.filterlistSE_2) <- str_replace(names(lb.filterlistSE_2),".*salmon/","")
  
  part1 <- suppressMessages(round(tximport(files = lb.filterlistSE_1, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  part2 <- suppressMessages(round(tximport(files = lb.filterlistSE_2, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  part3 <- suppressMessages(round(tximport(files = lb.filterlistSE_3, 
                                           type = "salmon",txOut=TRUE)$counts))
  
  all(colnames(part1)==colnames(part2))
  countsSE <- part1 + part2
  countsSE <- cbind(part3, countsSE)
  
  ###complete SE counts and sample processing
  ###merge counts and samples of ZE and SE into new set of sample and count files
  
  samples_SE_plus_ZE <- bind_rows(samplesZE_Final,samplesSE)
  samples_SE_plus_ZE$Tissue <- factor(samples_SE_plus_ZE$Tissue)
  samples_SE_plus_ZE$Experiment <- factor(samples_SE_plus_ZE$Experiment)
  samples_SE_plus_ZE$Time <- factor(samples_SE_plus_ZE$Time)
  
  counts_SE_plus_ZE <- cbind(countsZE_Final,countsSE)
  }
  
  #TT Model Modification
  {
    #Unique model using only B4-B6 time points, with SE and ZE tissue, using the Tissue * Time Model.
    #remove all samples with times that are not B4-B6, and remove FMG tissue
    samples.sz.ttmodel <- subset(samples_SE_plus_ZE, samples_SE_plus_ZE$Tissue != "FMG")
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B0"))
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B1"))
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B2"))
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B3"))
    #samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B7"))
    #samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B8"))
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B9"))
    samples.sz.ttmodel <- subset(samples.sz.ttmodel, samples.sz.ttmodel$Time != c("B10"))
    #do the removal after the model... Keep only times that are seen in both SE and ZE
    
    samples.sz.ttmodel <- subset(samples.sz.ttmodel,samples.sz.ttmodel$NGI.ID != "P464_205")
    counts.sz.ttmodel <- counts_SE_plus_ZE
    counts.sz.ttmodel <- counts.sz.ttmodel[ ,samples.sz.ttmodel$NGI.ID]
    colnames(counts.sz.ttmodel) == samples.sz.ttmodel$NGI.ID
  }
}

#SE_Germ Dataset
{
  lb.filelistSE_Germ <- list.files("/mnt/picea/projects/spruce/uegertsdotter/SE-germinants/salmon", 
                                          recursive = TRUE, 
                                          pattern = "quant.sf",
                                          full.names = TRUE)
  
  #################start of SomaticEmbGerm
  ##splitting data in order to sum the counts
  lb.filterlistSE_Germ <- lb.filelistSE_Germ
  lb.filterlistSE_Germ_1 <- (str_subset(lb.filterlistSE_Germ, "AC7", negate = TRUE))
  names(lb.filterlistSE_Germ_1) <- sapply(lapply(strsplit(lb.filterlistSE_Germ_1,"_"),"[",4:5),paste,collapse="_")
  lb.filterlistSE_Germ_2 <- (str_subset(lb.filterlistSE_Germ, "BC7", negate = TRUE))
  names(lb.filterlistSE_Germ_2) <- sapply(lapply(strsplit(lb.filterlistSE_Germ_2,"_"),"[",4:5),paste,collapse="_")
  
  
  part1 <- suppressMessages(round(tximport(files = lb.filterlistSE_Germ_1, 
                                           type = "salmon",txOut=TRUE)$counts))
  part2 <- suppressMessages(round(tximport(files = lb.filterlistSE_Germ_2, 
                                           type = "salmon",txOut=TRUE)$counts))
  all(colnames(part1)==colnames(part2))
  countsSE_Germ <- part1 + part2
  samplesSE_Germ
  #################end of EmbGerm
}

#ZE Only
samplesZE_Final
lb.g_ZE <- countsZE_Final
#export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(lb.g_ZE,file=here("analysis/salmon/ZE-ZF-unnormalised-gene-expression_data.csv"))

#4Datasets
samples_4sets <- bind_rows(samplesZE,samples29SEED,samplesSE_Germ,samplesSE)
lb.g_4sets <- cbind(countsZE, counts29Seed, countsSE_Germ, countsSE)
#export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(lb.g_4sets,file=here("analysis/salmon/4Datasets-unnormalised-gene-expression_data.csv"))

#Combine ZE and 29Seed for Batch Correction
samples_29z <- bind_rows(samplesZE,samples29SEED)
lb.g_29z <- cbind(countsZE, counts29Seed)
#export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(lb.g_29z,file=here("analysis/salmon/ZE-29Seed-unnormalised-gene-expression_data.csv"))

#Combine SE adn ZE for later use in Batch correction application
#SE vs ZE
samples_SE_plus_ZE
samples_sz <- bind_rows(samplesZE, samplesSE)
samples_sz <- samples_SE_plus_ZE
lb.g_sz <- counts_SE_plus_ZE
#export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(lb.g_sz,file=here("analysis/salmon/ZE-SE-unnormalised-gene-expression_data.csv"))


#DDS library
#ZE Only
dds.ze <- DESeqDataSetFromMatrix(
  countData = lb.g_ZE,
  colData = samplesZE_Final,
  design = ~ Tissue + Time + Tissue:Time)
save(dds.ze,file=here("analysis/salmon/ZE-ZF-Dataset-dds.rda"))

#ZE vs 29Seed
dds.29z <- DESeqDataSetFromMatrix(
  countData = lb.g_29z,
  colData = samples_29z,
  design = ~ Experiment)
save(dds.29z,file=here("analysis/salmon/ZE-29Seed-dds.rda"))

#4Datasets
dds.4sets <- DESeqDataSetFromMatrix(
  countData = lb.g_4sets,
  colData = samples_4sets,
  design = ~ Experiment)
save(dds.4sets,file=here("analysis/salmon/4Datasets-dds.rda"))

lb.g_sz
samples_sz

#SE vs ZE
dds.sz <- DESeqDataSetFromMatrix(
  countData = lb.g_sz,
  colData = samples_sz,
  design = ~ Experiment)
save(dds.sz,file=here("analysis/salmon/ZE-SE-Dataset-dds.rda"))
save(dds.sz,file=here("analysis/salmon/ZE-SE-dds.rda"))

dds.sz.ttmodel <- DESeqDataSetFromMatrix(
  countData = counts.sz.ttmodel,
  colData = samples.sz.ttmodel,
  design = ~Tissue * Time)
#dds.sz.ttmodel <- dds.sz.ttmodel[,!(dds.sz.ttmodel$NGI.ID == "P464_205")]
#view(samples_SE_plus_ZE)
#subset(samples_SE_plus_ZE,samples_SE_plus_ZE$Experiment == "Somatic Embryogenesis")
save(dds.sz.ttmodel,file=here("analysis/salmon/ZvsS-B4-B6-TTModel-dds.rda"))
#dds.sz.ttmodel.de$NGI.ID

dds.ze
dds.29z
dds.4sets
dds.sz
dds.sz.ttmodel

#' #ZE DATASET ONLY
{
dds.ze <- estimateSizeFactors(dds.ze)
sizes <- sizeFactors(dds.ze)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
vsd.ze <- varianceStabilizingTransformation(dds.ze, blind=TRUE)
vst.ze <- assay(vsd.ze)
vst.ze <- vst.ze - min(vst.ze)
pc.ze <- prcomp(t(vst.ze))
pc.dat.ze <- bind_cols(PC1=pc.ze$x[,1],
                    PC2=pc.ze$x[,2],
                    samplesZE)
ggplot(pc.dat.ze,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))
interplot <- ggplot(pc.dat.ze,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds.ze <- factor(paste(samplesZE$Tissue,samplesZE$Time))
sels.ze <- rangeFeatureSelect(counts=vst.ze,
                           conditions=conds.ze,
                           nrep=3)
vstCutoff_ZE <- 4+1

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples.
#' It appears that generally there is little difference between the samples across all genes
#'  - a small difference is noticable between ZE tissues and FMG-S tissues, however 
#' generally the expression levels are relatively balanced apart from the one sample 
#' in FMG B4. Somatic samples had a small section of more highly expressed genes 
#' compared to ZE and FMG.
heatmap.2(t(scale(t(vst.ze[sels.ze[[vstCutoff_ZE]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)
}

#' #29Seed vs ZE
{
dds.29z <- estimateSizeFactors(dds.29z)
sizes <- sizeFactors(dds.29z)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
vsd.29z <- varianceStabilizingTransformation(dds.29z, blind=TRUE)
vst.29z <- assay(vsd.29z)
vst.29z <- vst.29z - min(vst.29z)
pc.29z <- prcomp(t(vst.29z))
pc.dat.29z <- bind_cols(PC1=pc.29z$x[,1],
                    PC2=pc.29z$x[,2],
                    samples_29z)
ggplot(pc.dat.29z,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))
interplot <- ggplot(pc.dat.29z,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds.29z <- factor(paste(samples_29z$Tissue,samples_29z$Time))
sels.29z <- rangeFeatureSelect(counts=vst.29z,
                           conditions=conds.29z,
                           nrep=3)
vstCutoff_29Z <- 4+1

#space for heatmap

#end heatmap

#batch correction
{
  mat <- assay(vsd.29z)
  x=mat
  colnames(x) <- dds$NGI.ID
  x <- x[,-(48:64), drop = FALSE]
  x <- x[,-(1:21), drop = FALSE]
  ###clip off from x, column 21 and below
  ###then, clip off 48 to 65. Keep the rest.
  batch = vsd.29z$Experiment
  batch <- batch[-(48:64)]
  batch <- batch[-(1:21)]
  length(batch)
  design = matrix(1, ncol(x), 1) 
  batch <- as.factor(batch)
  contrasts(batch) <- contr.sum(levels(batch))
  batch <- model.matrix(~batch)[, -1, drop = FALSE]
  X.batch <- batch
  fit <- lmFit(x, cbind(design, X.batch))
  #beta is the actual batch correction for each gene
  beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
  beta[is.na(beta)] <- 0
  
  ####use the above, once corrected using B8-B9 AND talking with Ioana or Ulrika about our suspicion that 29Seed were mature seeds
  # Use the computed batch effect to normalise ZE and SE
  # Change x to be vst (ZE + SE) and X.batch, extend the -1 to as many as SE
  ###here we use the model on ZE and 29Seed
  as.matrix(x) - beta %*% t(X.batch)
  fullbatch <- vsd.29z$Experiment
  fullbatch <- as.factor(fullbatch)
  contrasts(fullbatch) <- contr.sum(levels(fullbatch))
  fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]
  X.fullbatch <- fullbatch
  as.matrix(mat) - beta %*% t(X.fullbatch)
  
  
  vsd.29z <- as.matrix(mat) - beta %*% t(X.fullbatch)
  vst.29z <- vsd.29z
  vst.29z <- vst.29z - min(vst.29z)
}
#batch correction should be completed - redo pca to look at the difference

#end pca

}
#NOTICE - PERFORM 29SEED_ZYGOTIC CODE BEFORE PERFORMING SZ

#' #4SETS
{
dds.4sets <- estimateSizeFactors(dds.4sets)
sizes <- sizeFactors(dds.4sets)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
vsd.4sets <- varianceStabilizingTransformation(dds.4sets, blind=TRUE)
vst.4sets <- assay(vsd.4sets)
vst.4sets <- vst.4sets - min(vst.4sets)
pc.4sets <- prcomp(t(vst.4sets))
pc.dat.4sets <- bind_cols(PC1=pc.4sets$x[,1],
                    PC2=pc.4sets$x[,2],
                    samples_4sets)
ggplot(pc.dat.4sets,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))
interplot <- ggplot(pc.dat.4sets,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds.4sets <- factor(paste(samples_4sets$Tissue,samples_4sets$Time))
sels.4sets <- rangeFeatureSelect(counts=vst.4sets,
                           conditions=conds.4sets,
                           nrep=3)
vstCutoff_4sets <- 4+1

heatmap.2(t(scale(t(vst.4sets[sels.4sets[[vstCutoff_4sets]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)



}

#' #SZ EXPERIMENT MODEL
{
dds.sz <- estimateSizeFactors(dds.sz)
sizes <- sizeFactors(dds.sz)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
vsd.sz <- varianceStabilizingTransformation(dds.sz, blind=TRUE)
vst.sz <- assay(vsd.sz)
vst.sz <- vst.sz - min(vst.sz)
pc.sz <- prcomp(t(vst.sz))
pc.dat.sz <- bind_cols(PC1=pc.sz$x[,1],
                    PC2=pc.sz$x[,2],
                    samples_sz)
ggplot(pc.dat.sz,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))
interplot <- ggplot(pc.dat.sz,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds.sz <- factor(paste(samples_sz$Tissue,samples_sz$Time))
sels.sz <- rangeFeatureSelect(counts=vst.sz,
                           conditions=conds.sz,
                           nrep=3)
vstCutoff_SZ <- 4+1
}
#space for heatmap

#end heatmap

#batch correction - relies on 29z batch correction first

  mat_SE_ZE <- assay(vsd.sz)
  x_SE_ZE <- mat_SE_ZE
  
  length(colnames(vsd.sz))
  vsd.sz$Experiment
  
  as.matrix(x_SE_ZE) - beta %*% t(X.batch)
  fullbatch <- vsd.sz$Experiment
  fullbatch <- as.factor(fullbatch)
  contrasts(fullbatch) <- contr.sum(levels(fullbatch))
  fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]
  X.fullbatch <- fullbatch
  as.matrix(mat_SE_ZE) - beta %*% t(X.fullbatch)
  samples_sz$Tissue
  
  vsd.sz <- as.matrix(mat_SE_ZE) - beta %*% t(X.fullbatch)
  vst.sz <- vsd.sz
  vst.sz <- vst.sz - min(vst.sz)


#' #SZ TT MODEL
{
dds.sz.ttmodel <- estimateSizeFactors(dds.sz.ttmodel)
sizes <- sizeFactors(dds.sz.ttmodel)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
vsd.sz.ttmodel <- varianceStabilizingTransformation(dds.sz.ttmodel, blind=TRUE)
vst.sz.ttmodel <- assay(vsd.sz.ttmodel)
vst.sz.ttmodel <- vst.sz.ttmodel - min(vst.sz.ttmodel)
pc.sz.ttmodel <- prcomp(t(vst.sz.ttmodel))
pc.dat.sz.ttmodel <- bind_cols(PC1=pc.sz.ttmodel$x[,1],
                    PC2=pc.sz.ttmodel$x[,2],
                    samples.sz.ttmodel)
ggplot(pc.dat.sz.ttmodel,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))
interplot <- ggplot(pc.dat.sz.ttmodel,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds.sz.ttmodel <- factor(paste(samples.sz.ttmodel$Tissue,samples.sz.ttmodel$Time))
sels.sz.ttmodel <- rangeFeatureSelect(counts=vst.sz.ttmodel,
                           conditions=conds.sz.ttmodel,
                           nrep=3)
vstCutoff_SZ.ttmodel <- 4+1

#space for heatmap

#end heatmap

#batch correction - relies on 29z batch correction first

mat.sz.ttmodel <- assay(vsd.sz.ttmodel)
x.sz.ttmodel <- mat.sz.ttmodel

length(colnames(vsd.sz.ttmodel))
vsd.sz.ttmodel$Experiment

as.matrix(x.sz.ttmodel) - beta %*% t(X.batch)
fullbatch <- vsd.sz.ttmodel$Experiment
fullbatch <- as.factor(fullbatch)
contrasts(fullbatch) <- contr.sum(levels(fullbatch))
fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]
X.fullbatch <- fullbatch
as.matrix(mat.sz.ttmodel) - beta %*% t(X.fullbatch)
samples_sz$Tissue

vsd.sz.ttmodel <- as.matrix(mat.sz.ttmodel) - beta %*% t(X.fullbatch)
vst.sz.ttmodel <- vsd.sz.ttmodel
vst.sz.ttmodel <- vst.sz.ttmodel - min(vst.sz.ttmodel)

}













#Tissue Specificity Section


#Gopher Enrichment
  {
  metaDF <- data.frame(dds_SE_ZE$NGI.ID,dds_SE_ZE$Time,dds_SE_ZE$Tissue,dds_SE_ZE$Experiment)
  metaDF$NGI.ID <- as.character(metaDF$NGI.ID)
  colnames(metaDF) <- c("NGI.ID","Time","Tissue","Experiment")
  D3Time <- plot3Dpca(t(vst_SE_ZE_featureselected), metaDF, c("Time"),colors = cust_colors, title="ZE vs SE batch corrected (Time)", inverse=c(FALSE,FALSE,FALSE))
  D3Tissue <- plot3Dpca(t(vst_SE_ZE_featureselected), metaDF, c("Tissue"), title="ZE vs SE batch corrected (Tissue)", inverse=c(FALSE,FALSE,FALSE))
  }

got_pc <- get_pca(pc, element = c("var", "ind"))
top500 <- fviz_contrib(pc, choice = "var", axes = 1, top = 500)
top50 <- fviz_contrib(pc, choice = "var", axes = 1, top = 50)
toptop <- fviz_contrib(pc, choice = "var", axes = 1, top = 13346)
t<-top500
# see if you can extract the "red-line"

D1 <- got_pc$contrib[,1]
got_pc$cor
samples_SE_plus_ZE


# if not play with the contributions
# which(cumsum(sort(t$data$contrib,decreasing = TRUE))>=50)[1]
# 
(t$data$contrib - mean(t$data$contrib)) > 0
which(cumsum(sort(t$data$contrib,decreasing = TRUE))>=50)[1]
which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1] #485 genes
which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=100)[1] #13346 genes

##this gives me the sorted genes in decreasing order - know that the first 485 give me 50%, only keep the ones above 485
genes <- sort(got_pc$contrib[,1],decreasing = TRUE)
#genes that give 50% of PC1
which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1] #485 genes
which(cumsum(sort(got_pc$contrib[,2],decreasing = TRUE))>=50)[1] #526 genes
which(cumsum(sort(got_pc$contrib[,3],decreasing = TRUE))>=50)[1] #592 genes
genes_PC1_50 <- genes[1:485]
#genes that give 10% of PC1
which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=10)[1] #485 genes
genes_PC1_10 <- genes[1:which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=10)[1]]

genes_PC1_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1]]
genes_PC2_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,2],decreasing = TRUE))>=50)[1]]
genes_PC3_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,3],decreasing = TRUE))>=50)[1]]

genes_PC1_10[1]
names(genes_PC1_10[1])
# look at line plots
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

line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[1]))
line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[2]))
line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[4]))
line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[12]))

genes_PC1_10_names <- names(genes_PC1_10)
genes_PC1_50_names <- names(genes_PC1_50)
genes_PC2_50_names <- names(genes_PC2_50)
genes_PC3_50_names <- names(genes_PC3_50)

background <- rownames(vst_SE_ZE_featureselected)
background <- str_replace(background,"[.]1","")
genes_PC1_10_names <- str_replace(genes_PC1_10_names,"[.]1","")
rownames(vst_SE_ZE_featureselected)
enr <- gopher(genes = genes_PC1_10_names, background = background,task = "go",url = "pabies")

#processing gene_names so that they remove the .1 at the end of every gene name identity (gopher cannot use it if its there)
genes_PC1_50_names <- names(genes_PC1_50)
genes_PC1_50_names <- str_replace(genes_PC1_50_names,"[.]1","")
genes_PC2_50_names <- str_replace(genes_PC2_50_names,"[.]1","")
genes_PC3_50_names <- str_replace(genes_PC3_50_names,"[.]1","")

#normal gopher, with featureselected background
enr <- gopher(genes = genes_PC1_50_names, background = background,task = c("go","mapman"),url = "pabies")
enr2 <- gopher(genes = genes_PC2_50_names, background = background,task = c("go","mapman"),url = "pabies")
enr3 <- gopher(genes = genes_PC3_50_names, background = background,task = c("go","mapman"),url = "pabies")

##testing without the featureselected background
enrnoback <- gopher(genes = genes_PC1_50_names, background = NULL,task = c("go","mapman"),url = "pabies")
enr2noback <- gopher(genes = genes_PC2_50_names, background = NULL,task = c("go","mapman"),url = "pabies")
enr3noback <- gopher(genes = genes_PC3_50_names, background = NULL,task = c("go","mapman"),url = "pabies")

# gopher
source(here("UPSCb-common/src/R/gopher.R"))
enr <- gopher(genes = "contrib", background = "the genenames selected fro featureSelect",task = c("go","mapman"),url = "pabies")

enr2$go$id
enr3$go$id

commons1 <- str_subset(enr$go$id,enr2$go$id)
commons2 <- str_subset(enr$go$id,enr3$go$id)
commons3 <- str_subset(enr2$go$id,enr3$go$id)

commons1
commons2
commons3
str_subset(commons1, commons3)
str_subset(commons1, commons2)
str_subset(commons2, commons3)

# extract the ID,padj columns and perform a visualization in REVIGO (or ask Alonso for his code to plot treemaps)
enr_exfil <- subset(enr$go, select = c(id, padj))
genes_PC1_50_names
enr_exfil
bind_cols(unlist(genes_PC1_50_names),enr_exfil)

source(here("Rtoolbox/src/plotEnrichedTreemap.R"))

genes_PC1_50_names
enr$go$id
enr2
enr3

#PC1 GO and Mapman Treemap
plotEnrichedTreemap(enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr, enrichment = "mapman", clusterColor = "#1B75BC")

#PC2 GO and Mapman Treemap
plotEnrichedTreemap(enr2, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr2, enrichment = "mapman", clusterColor = "#1B75BC")

#PC3 GO and Mapman Treemap
plotEnrichedTreemap(enr3, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr3, enrichment = "mapman", clusterColor = "#1B75BC")

#PC1 GO and Mapman Noback
plotEnrichedTreemap(enrnoback, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enrnoback, enrichment = "mapman", clusterColor = "#1B75BC")

#PC2 GO and Mapman Noback
plotEnrichedTreemap(enr2noback, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr2noback, enrichment = "mapman", clusterColor = "#1B75BC")

#PC3 GO and Mapman Noback
plotEnrichedTreemap(enr3noback, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr3noback, enrichment = "mapman", clusterColor = "#1B75BC")





enr$go$name
enr2$go$name
enr3$go$name

enrnoback$go$name
enr2noback$go$name
enr3noback$go$name

str_subset(enrnoback$go$name, enr3noback$go$name)



samples_SE_plus_ZE$NGI.ID

####redo analysis, but remove FMG completely from the equation.
{
  samples_SE_plus_ZE_minusFMG <- subset(samples_SE_plus_ZE, samples_SE_plus_ZE$Tissue != "FMG")
  ##delete from counts and samples_SE_plus_ZE
  dds_SE_ZE_minusFMG <- dds_SE_ZE_minusFMG[,(dds_SE_ZE_minusFMG$Tissue == "FMG")]
  
  
  dds_SE_ZE_minusFMG <- dds_SE_ZE_minusFMG[,(dds_SE_ZE_minusFMG$NGI.ID == samples_SE_plus_ZE_minusFMG$NGI.ID)]
  dds_SE_ZE_minusFMG$NGI.ID
  
  vsd_SE_ZE_minusFMG <- varianceStabilizingTransformation(dds_SE_ZE_minusFMG,blind=FALSE)
  vst_SE_ZE_minusFMG <- assay(vsd_SE_ZE_minusFMG)
  vst_SE_ZE_minusFMG <- vst_SE_ZE_minusFMG - min(vst_SE_ZE_minusFMG)
  #' Run the feature selection
  #' Filter for noise
  #' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
  conds <- factor(paste(dds_SE_ZE_minusFMG$Tissue,dds_SE_ZE_minusFMG$Time))
  sels <- rangeFeatureSelect(counts=vst_SE_ZE_minusFMG,
                             conditions=conds,
                             nrep=3)
  
  vstCutoff <- 4+1
  vst_SE_ZE_minusFMG #66056 rows cutoff
  vst_SE_ZE_minusFMG[sels[[vstCutoff]],] #42544 rows cutoff
  
  vst_SE_ZE_minusFMG_featureselected <- vst_SE_ZE_minusFMG[sels[[vstCutoff]],]
  
  pc <- prcomp(t(vst_SE_ZE_minusFMG_featureselected))
  #' ### 2D
  {
    pc.dat <- bind_cols(PC1=pc$x[,1],
                        PC2=pc$x[,2],
                        samples_SE_plus_ZE_minusFMG)
    
    ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue)) + 
      geom_point(size=2) + 
      ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
      scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
      scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))
    
    #' ### Interactive PCA Plot
    suppressPackageStartupMessages(library(plotly))
    
    interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment,text=NGI.ID)) +
      geom_point(size=2) +
      ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
    
    dds_SE_ZE_minusFMG$Tissue
    
    
    ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                    yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
    
    metaDF_minusFMG <- data.frame(dds_SE_ZE_minusFMG$NGI.ID,dds_SE_ZE_minusFMG$Time,dds_SE_ZE_minusFMG$Tissue,dds_SE_ZE_minusFMG$Experiment)
    colnames(metaDF_minusFMG) <- c("NGI.ID","Time","Tissue","Experiment")
    metaDF_minusFMG$dds_SE_ZE_minusFMG.NGI.ID <- as.character(metaDF_minusFMG$NGI.ID)
    
    D3Time <- plot3Dpca(t(vst_SE_ZE_minusFMG_featureselected), metaDF_minusFMG, c("Time"),colors = cust_colors, title="ZE vs SE batch corrected (Time)", inverse=c(FALSE,FALSE,FALSE))
    D3Tissue <- plot3Dpca(t(vst_SE_ZE_minusFMG_featureselected), metaDF_minusFMG, c("Tissue"), title="ZE vs SE batch corrected (Tissue)", inverse=c(FALSE,FALSE,FALSE))
    
    D3Time
    D3Tissue
  }
}




{
  got_pc <- get_pca(pc, element = c("var", "ind"))
  
  # if not play with the contributions
  # which(cumsum(sort(t$data$contrib,decreasing = TRUE))>=50)[1]
  # 
  
  ##this gives me the sorted genes in decreasing order - know that the first 485 give me 50%, only keep the ones above 485
  genes <- sort(got_pc$contrib[,1],decreasing = TRUE)
  #genes that give 50% of PC1
  which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1] #485 genes
  which(cumsum(sort(got_pc$contrib[,2],decreasing = TRUE))>=50)[1] #526 genes
  which(cumsum(sort(got_pc$contrib[,3],decreasing = TRUE))>=50)[1] #592 genes
  genes_PC1_50 <- genes[1:485]
  #genes that give 10% of PC1
  which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=10)[1] #485 genes
  genes_PC1_10 <- genes[1:which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=10)[1]]
  
  genes <- sort(got_pc$contrib[,1],decreasing = TRUE)
  genes_PC1_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,1],decreasing = TRUE))>=50)[1]]
  genes <- sort(got_pc$contrib[,2],decreasing = TRUE)
  genes_PC2_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,2],decreasing = TRUE))>=50)[1]]
  genes <- sort(got_pc$contrib[,3],decreasing = TRUE)
  genes_PC3_50 <- genes[1:which(cumsum(sort(got_pc$contrib[,3],decreasing = TRUE))>=50)[1]]
  
  length(genes_PC1_50)
  length(genes_PC2_50)
  length(genes_PC3_50)
  
  genes_PC1_10[1]
  names(genes_PC1_10[1])
  # look at line plots
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
  
  line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[1]))
  line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[2]))
  line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[4]))
  line_plot(dds_SE_ZE,vst_SE_ZE,names(genes_PC1_10[12]))
  
  genes_PC1_50_names <- names(genes_PC1_50)
  genes_PC2_50_names <- names(genes_PC2_50)
  genes_PC3_50_names <- names(genes_PC3_50)
  
  background <- rownames(vst_SE_ZE_minusFMG_featureselected)
  background <- str_replace(background,"[.]1","")
  genes_PC1_10_names <- str_replace(genes_PC1_10_names,"[.]1","")
  rownames(vst_SE_ZE_minusFMG_featureselected)
  enr <- gopher(genes = genes_PC1_10_names, background = background,task = "go",url = "pabies")
  
  #processing gene_names so that they remove the .1 at the end of every gene name identity (gopher cannot use it if its there)
  genes_PC1_50_names <- names(genes_PC1_50)
  genes_PC1_50_names <- str_replace(genes_PC1_50_names,"[.]1","")
  genes_PC2_50_names <- str_replace(genes_PC2_50_names,"[.]1","")
  genes_PC3_50_names <- str_replace(genes_PC3_50_names,"[.]1","")
  
  # gopher
  source(here("UPSCb-common/src/R/gopher.R"))
  enr <- gopher(genes = "contrib", background = "the genenames selected fro featureSelect",task = c("go","mapman"),url = "pabies")
  
  #normal gopher, with featureselected background
  enr <- gopher(genes = genes_PC1_50_names, background = background,task = c("go","mapman"),url = "pabies")
  enr2 <- gopher(genes = genes_PC2_50_names, background = background,task = c("go","mapman"),url = "pabies")
  enr3 <- gopher(genes = genes_PC3_50_names, background = background,task = c("go","mapman"),url = "pabies")
  
  ##testing without the featureselected background
  enrnoback <- gopher(genes = genes_PC1_50_names, background = NULL,task = c("go","mapman"),url = "pabies")
  enr2noback <- gopher(genes = genes_PC2_50_names, background = NULL,task = c("go","mapman"),url = "pabies")
  enr3noback <- gopher(genes = genes_PC3_50_names, background = NULL,task = c("go","mapman"),url = "pabies")
  
  
  
  enr2$go$id
  enr3$go$id
  
  commons1 <- str_subset(enr$go$id,enr2$go$id)
  commons2 <- str_subset(enr$go$id,enr3$go$id)
  commons3 <- str_subset(enr2$go$id,enr3$go$id)
  
  commons1
  commons2
  commons3
  str_subset(commons1, commons3)
  str_subset(commons1, commons2)
  str_subset(commons2, commons3)
  
  # extract the ID,padj columns and perform a visualization in REVIGO (or ask Alonso for his code to plot treemaps)
  enr_exfil <- subset(enr$go, select = c(id, padj))
  genes_PC1_50_names
  enr_exfil
  bind_cols(unlist(genes_PC1_50_names),enr_exfil)
  
  source(here("Rtoolbox/src/plotEnrichedTreemap.R"))
  
  genes_PC1_50_names
  enr$go$id
  enr2
  enr3
  
  #PC1 GO and Mapman Treemap
  plotEnrichedTreemap(enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #PC2 GO and Mapman Treemap
  plotEnrichedTreemap(enr2, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enr2, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #PC3 GO and Mapman Treemap
  plotEnrichedTreemap(enr3, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enr3, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #PC1 GO and Mapman Noback
  plotEnrichedTreemap(enrnoback, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enrnoback, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #PC2 GO and Mapman Noback
  plotEnrichedTreemap(enr2noback, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enr2noback, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #PC3 GO and Mapman Noback
  plotEnrichedTreemap(enr3noback, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(enr3noback, enrichment = "mapman", clusterColor = "#1B75BC")
  
  
  
  
  
  enr$go$name
  enr2$go$name
  enr3$go$name
  
  enr$mapman$name
  enr2$mapman$name
  enr3$mapman$name
  
  enrnoback$go$name
  enr2noback$go$name
  enr3noback$go$name
  
}



