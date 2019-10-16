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

### Sample preprocessing
##need to delete User.ID, Sample.ID, Replicate, Mreads and X..Q30 (columns 2,4,6,7,8)
samples <- filter(samples, !grepl("P11562_112",NGI.ID))
samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )
samples
samplesZE <- subset(samples, samples$Experiment == "Zygotic Embryogenesis")
samplesSEED <- subset(samples, samples$Experiment == "29Seed")
samplesSE <- subset(samples, samples$Experiment == "Somatic Embryogenesis")




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

lb.g <- cbind(countsZE_Final, counts29Seed)
stopifnot(all(str_which(colnames(lb.g), samples$NGI.ID) == 1:length(colnames(lb.g))))

counts <- lb.g

################################SE IMPORTING
#' * Metadata
#' Sample information ########### need sample info?
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

###complete SE counts and sample processing
###merge counts and samples of ZE and SE into new set of sample and count files

samples_SE_plus_ZE <- bind_rows(samplesZE_Final,samplesSE)
samples_SE_plus_ZE$Tissue <- factor(samples_SE_plus_ZE$Tissue)
samples_SE_plus_ZE$Experiment <- factor(samples_SE_plus_ZE$Experiment)
samples_SE_plus_ZE$Time <- factor(samples_SE_plus_ZE$Time)

counts_SE_plus_ZE <- cbind(countsZE_Final,countsSomaticEmb)

length(colnames(counts_SE_plus_ZE))
length(samples_SE_plus_ZE$NGI.ID)

stopifnot(all(str_which(colnames(counts_SE_plus_ZE), samples_SE_plus_ZE$NGI.ID) == 1:length(colnames(counts_SE_plus_ZE))))


####ZE and 29SEED FIRST
#' ## Raw Data QC analysis
#' Check how many genes are never expressed - reasonable level of non-expressed genes indicated.
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
       aes(x,y,col=Experiment,fill=Tissue)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())


ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(lb.filterlist),samples$NGI.ID),]),
       aes(x,y,col=Experiment,fill=)) + geom_col() + 
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
  geom_density()  + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Tissue=samples$Tissue[match(ind,samples$NGI.ID)]) %>% 
  mutate(Experiment=samples$Experiment[match(ind,samples$NGI.ID)])

#' Color by Experiment
ggplot(dat,aes(x=values,group=ind,col=Tissue)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Time
ggplot(dat,aes(x=values,group=ind,col=Experiment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("analysis/salmon/ZE-29Seed-unnormalised-gene-expression_data.csv"))
############## change export location name


#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#' 
#' This first dds is the ZE and 29Seed dds file, following that, the next dds file will be for SE and ZE
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~Experiment)
samples$Time
dds <- dds[,!(dds$NGI.ID == "P11562_148")]
colnames(dds) <- dds$NGI.ID
dds$NGI.ID
samples <- subset(samples, samples$NGI.ID != "P11562_148")
save(dds,file=here("analysis/salmon/ZE-29Seed-dds.rda"))

length(colnames(counts))
length(samples$NGI.ID)

###adjusting processing
#remove S (B1-B3) and the B10 sample
backup_SE_ZE <- samples_SE_plus_ZE
#samples_SE_plus_ZE <- backup_SE_ZE
samples_SE_plus_ZE$Time
destroy <- samples_SE_plus_ZE
destroy <- subset(samples_SE_plus_ZE,samples_SE_plus_ZE$Experiment == "Zygotic Embryogenesis")
destroy_SE <- subset(samples_SE_plus_ZE,samples_SE_plus_ZE$Experiment == "Somatic Embryogenesis")

destroy <- subset(destroy,destroy$Time != "B1")
destroy <- subset(destroy,destroy$Time != "B2")
destroy <- subset(destroy,destroy$Time != "B3")
destroy <- subset(destroy,destroy$Time != "B10")
destroy$Time

samples_SE_plus_ZE <- bind_rows(destroy, destroy_SE)
##delete from counts and samples_SE_plus_ZE
counts_SE_plus_ZE <- counts_SE_plus_ZE[ ,samples_SE_plus_ZE$NGI.ID]


dds_SE_ZE <- DESeqDataSetFromMatrix(
  countData = counts_SE_plus_ZE,
  colData = samples_SE_plus_ZE,
  design = ~Experiment)
samples_SE_plus_ZE$Time
dds_SE_ZE <- dds_SE_ZE[,!(dds_SE_ZE$NGI.ID == "P11562_148")]
colnames(dds_SE_ZE) <- dds_SE_ZE$NGI.ID
dds$NGI.ID
samples_SE_plus_ZE <- subset(samples_SE_plus_ZE,samples_SE_plus_ZE$NGI.ID != "P11562_148")
save(dds_SE_ZE,file=here("analysis/salmon/ZE-SE-dds.rda"))

load(here("analysis/salmon/ZE-29Seed-dds.rda"))
load(here("analysis/salmon/ZE-SE-dds.rda"))

suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(bladderbatch))
suppressPackageStartupMessages(library(pamr))
suppressPackageStartupMessages(library(limma))

library(convert)
library(Biobase)
require(Biobase)
source(here("src/R/ComBat-seq-master/ComBat_seq.R"))
source(here("src/R/ComBat-seq-master/helper_seq.R"))
library(edgeR)
#' ## Normalisation for visualisation
#' the normalisation is aware to take advantage of the model to determine the dispersion
#' 
#' ZE and 29Seed First, SE and ZE after
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

vsd_SE_ZE <- varianceStabilizingTransformation(dds_SE_ZE,blind=FALSE)
vst_SE_ZE <- assay(vsd_SE_ZE)
vst_SE_ZE <- vst_SE_ZE - min(vst_SE_ZE)

#limma
source(here("UPSCb-common/src/R/percentile.R"))

#using ZE and 29Seed first, to get the batch effect
###mad is vsd_adjusted
mat <- assay(vsd)

###original full model - not the one to use - see below
{
x=mat
batch = vsd_adjusted$Experiment
design = matrix(1, ncol(x), 1) 
batch <- as.factor(batch)
contrasts(batch) <- contr.sum(levels(batch))
batch <- model.matrix(~batch)[, -1, drop = FALSE]
X.batch <- batch
fit <- lmFit(x, cbind(design, X.batch))
beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
beta[is.na(beta)] <- 0
samples$Time
}
# TODO ask Ulrika / Ioana Gaboreanu
# Redo the batch correction using only B8-B9 (guess is that it will be more accurate)
###trim the samples to remove the samples that are below B8, maybe remove above B10.
{
x=mat
colnames(x) <- dds$NGI.ID
x <- x[,-(48:64), drop = FALSE]
x <- x[,-(1:21), drop = FALSE]
###clip off from x, column 21 and below
###then, clip off 48 to 65. Keep the rest.
batch = vsd_adjusted$Experiment
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
fullbatch <- vsd$Experiment
fullbatch <- as.factor(fullbatch)
contrasts(fullbatch) <- contr.sum(levels(fullbatch))
fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]
X.fullbatch <- fullbatch
as.matrix(mat) - beta %*% t(X.fullbatch)


vsd <- as.matrix(mat) - beta %*% t(X.fullbatch)
vst <- vsd
vst <- vst - min(vst)
}
############################################################################################################################
#plot PCA ZE and 29Seed
{#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))

percent <- round(summary(pc)$importance[2,]*100)

samples$NGI.ID
samplesfinal <- samples_SE_plus_ZE

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samplesfinal)

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))

interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")




ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
}
#########################################################################################################################
###here we use the model on ZE and SE
{
mat_SE_ZE <- assay(vsd_SE_ZE)
x_SE_ZE <- mat_SE_ZE

length(colnames(vsd_SE_ZE))
vsd_SE_ZE$Experiment

as.matrix(x_SE_ZE) - beta %*% t(X.batch)
fullbatch <- vsd_SE_ZE$Experiment
fullbatch <- as.factor(fullbatch)
contrasts(fullbatch) <- contr.sum(levels(fullbatch))
fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]
X.fullbatch <- fullbatch
as.matrix(mat_SE_ZE) - beta %*% t(X.fullbatch)
samples_SE_plus_ZE$Tissue

vsd_SE_ZE <- as.matrix(mat_SE_ZE) - beta %*% t(X.fullbatch)
vst_SE_ZE <- vsd_SE_ZE
vst_SE_ZE <- vst_SE_ZE - min(vst_SE_ZE)

# plot PCA ZE vs. SE
pc <- prcomp(t(vst_SE_ZE))

percent <- round(summary(pc)$importance[2,]*100)

samples$NGI.ID
samplesfinal <- samples_SE_plus_ZE

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samplesfinal)

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))

interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")




ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))
}





















#####using LIMMA unaltered
dds_limma <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Experiment)
dds_limma$NGI.ID <- samples$NGI.ID

dds_limma <- dds_limma[,!(dds_limma$NGI.ID == "P11562_148")]
vsd_limma <- varianceStabilizingTransformation(dds_limma,blind=FALSE)
#############################
#############################
#############################
##attempting limma batch correction
mat <- assay(vsd_limma)
mat <- limma::removeBatchEffect(mat, vsd_limma$Experiment)
assay(vsd_limma) <- mat
vst_limma <- assay(vsd_limma)
vst_limma <- vst_limma - min(vst_limma)


colnames(vst_limma) <- colnames(dds_limma)

#####MERGE VST LIMMA AND VST HERE
vst <- vst_limma

meanSdPlot(vst[rowSums(counts_limma)>0,])
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))

percent <- round(summary(pc)$importance[2,]*100)
samples$NGI.ID
samplesfinal <- samples[!(samples$NGI.ID == "P11562_148"),]

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samplesfinal)

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Interactive PCA Plot
suppressPackageStartupMessages(library(plotly))

interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")




ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                                                yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

##########Limma looks like the best
##########apply this data, and put it together with SE.


















