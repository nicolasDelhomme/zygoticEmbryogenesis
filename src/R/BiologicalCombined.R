#' ---
#' title: "All SE and ZE datasets combined"
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
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(plotly)
  library(RColorBrewer)
  library(scatterplot3d)
  library(tidyverse)
  library(tximport)
  library(vsn)
  library(here)
})

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

##need to delete User.ID, Sample.ID, Replicate, Mreads and X..Q30 (columns 2,4,6,7,8)
samples <- subset(samples, select = -c(User.ID, Sample.ID, Replicate, Mreads, X..Q30) )
samples$Experiment

sZE <- subset(samples,samples$Experiment == "Zygotic Embryogenesis") 
s29 <- subset(samples,samples$Experiment == "29Seed")
sSG <- subset(samples,samples$Experiment == "Somatic Embryogenesis Germinants")
sSE <- subset(samples,samples$Experiment == "Somatic Embryogenesis")

sSE$Time <- sSE$Tissue
sSE$Time <- str_c("B",sSE$Time)
sSE$Tissue <- c("SE")

samples <- bind_rows(sZE,s29,sSG,sSE)
samples$Tissue <- factor(samples$Tissue)

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

#################end of ZE









lb.filelistSomaticEmbGerm <- list.files("/mnt/picea/projects/spruce/uegertsdotter/SE-germinants/salmon", 
                            recursive = TRUE, 
                            pattern = "quant.sf",
                            full.names = TRUE)

#################start of SomaticEmbGerm
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

#################end of EmbGerm



lb.filelist29Seed <- list.files("/mnt/picea/projects/spruce/uegertsdotter/29_Spruce_Seeds_Project/salmon", 
                            recursive = TRUE, 
                            pattern = "quant.sf",
                            full.names = TRUE)

lb.filterlist29Seed <- lb.filelist29Seed
names(lb.filterlist29Seed) <- sapply(lapply(strsplit(lb.filterlist29Seed,"_"),"[",7:8),paste,collapse="_")
counts29Seed <- suppressMessages(round(tximport(files = lb.filterlist29Seed, 
                                   type = "salmon",txOut=TRUE)$counts))
###need to recorder the ones with a letter infront of them


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
################# end of SomaticEmb



#' Filter and select only one duplicate from each sample.













#removed P11562_112 row from sample list
samples <- filter(samples, !grepl("P11562_112",NGI.ID))

lb.filterlist <- c(lb.filterlistZE,lb.filterlist29Seed,lb.filterlistSomaticEmbGerm1,lb.filterlistSomaticEmb)
lb.filterlist <- c(lb.filterlistZE,lb.filterlist29Seed,lb.filterlistSomaticEmbGerm1,lb.filterlistSomaticEmb)
#####Filelists have been combined, however not able to proceed - NGI IDs do not match.
#####All samples should have the same Tissue Type in the Somatic and 29Seed - Time is not changed either - they should all be the same.


stopifnot(all(str_which(lb.filterlist, samples$NGI.ID) == 1:length(lb.filterlist)))
#print(str_which(lb.filterlist, samples$NGI.ID) == 1:length(lb.filterlist))
##assign names to the filtered filelist (which removed L002)
names(lb.filterlist) <- samples$NGI.ID
#lb.filterlist <- lb.filterlist[samples$Tissue %in% c("ZE","FMG","S","Normal","Aberrant","Non-EMB","PEM","DKM","SM","LSM","SD","ED","RO","PLS","ZE-R")] ##Are we only looking at ZE?
#lb.filterlist <- lb.filterlist[samples$Experiment %in% c("Zygotic","29Seed","Somatic Embryogenesis Germinants","Somatic Embryogenesis")]

#' Read the expression at the gene level (there is one transcript per gene)
#lb.g <- suppressMessages(tximport(files = lb.filterlist, 
#                                  type = "salmon",txOut=TRUE))
#####redo this lb.g <- cbind(countsZE, countsSEED, countsSomaticEmbGerm, countsSomaticEmb) in the correct order
#####the order is ZE, Seeds, SomaticEmbGerm, SomaticEmb
lb.g <- cbind(countsZE, counts29Seed, countsEmbGerm, countsSomaticEmb)


counts <- lb.g













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
write.csv(counts,file=here("analysis/salmon/4Datasets-unnormalised-gene-expression_data.csv"))
############## change export location name


########################################want to get here to the Data normalisation
#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Experiment)

dds <- dds[,!(dds$NGI.ID == "P11562_148")]


save(dds,file=here("analysis/salmon/4Datasets-dds.rda"))

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

#' ### 3 first dimensions
#' Seems that different time points form small clusters and ZE and FMG tissue types appear to separate. These are Comp1 and Comp2 
#' which explains the different between most of the sampels except for one B4 ZE sample. Appears to be an outlier.
#' This seems to indicate that the Tissue and Time components explain the difference between samples.
#mar=c(5.1,4.1,4.1,2.1)
#scatterplot3d(pc$x[,1],
#              pc$x[,2],
#              pc$x[,3],
#              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
#              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
#              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
#              color=pal[as.integer(samples$Tissue)],
#              pch=c(17,19,15,13)[as.integer(samples$Experiment)])
#legend("topleft",
#       fill=pal[1:nlevels(samples$Tissue)],
#       legend=levels(samples$Tissue))
#legend("bottomright",
#       pch=c(17,19,15,13),
#       legend=levels(samples$Experiment))
#par(mar=mar)
samples$NGI.ID
samples <- samples[!(samples$NGI.ID == "P11562_148"),]

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Tissue,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Interactive PCA Plot
interplot <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Tissue,shape=Experiment,text=NGI.ID)) +
  geom_point(size=2) +
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")




ggplotly(interplot, tooltip = "all") %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                       yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))



#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 32000 genes - is this adequate for the QA?
conds <- factor(paste(samples$Tissue,samples$Experiment))
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vstCutoff <- 6+1

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples.
#' It appears that generally there is little difference between the samples across all genes
#'  - a small difference is noticable between ZE tissues and FMG-S tissues, however 
#' generally the expression levels are relatively balanced apart from the one sample 
#' in FMG B4. Somatic samples had a small section of more highly expressed genes 
#' compared to ZE and FMG.
hm <- heatmap.2(t(scale(t(vst[sels[[vstCutoff]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal,cexCol = .2)

plot(as.hclust(hm$colDendrogram))




#' Perform the VST, aware
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' # Feature selection
#' A cutoff at 1 VST seems appropriate for both
sels <- rangeFeatureSelect(
  counts=as.matrix(vst),
  conditions=colData(dds)$Experiment,nrep=3)

# check and decide on cutoff
cutoff=4
  sel <- sels[[cutoff]]

#' # Export
dir.create(here("data/seidr"),showWarnings=FALSE,recursive=TRUE)
           
           #' * gene by column, without names matrix
           write.table(t(vst[sel,]),
                       file=here("data/seidr/headless.tsv"),
                       col.names=FALSE,
                       row.names=FALSE,
                       sep="\t",quote=FALSE)
           
           #' * gene names, one row
           # TODO check IDs		   
           write.table(t(sub("\\.1$","",rownames(vst)[sel])),
                       file=here("data/seidr/genes.tsv"),
                       col.names=FALSE,
                       row.names=FALSE,
                       sep="\t",quote=FALSE)



#' ## Conclusion
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
