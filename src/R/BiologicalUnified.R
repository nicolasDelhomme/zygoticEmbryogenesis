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
samplesSEED <- subset(samples, samples$Experiment == "29Seed")
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
  countsSomaticEmb <- part1 + part2
  countsSomaticEmb <- cbind(part3, countsSomaticEmb)
  
  ###complete SE counts and sample processing
  ###merge counts and samples of ZE and SE into new set of sample and count files
  
  samples_SE_plus_ZE <- bind_rows(samplesZE_Final,samplesSE)
  samples_SE_plus_ZE$Tissue <- factor(samples_SE_plus_ZE$Tissue)
  samples_SE_plus_ZE$Experiment <- factor(samples_SE_plus_ZE$Experiment)
  samples_SE_plus_ZE$Time <- factor(samples_SE_plus_ZE$Time)
  
  counts_SE_plus_ZE <- cbind(countsZE_Final,countsSomaticEmb)
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
  countsEmbGerm <- part1 + part2
  samplesSE_Germ
  #################end of EmbGerm
}

#4Datasets
#Combine all the datasets


#Combine ZE and 29Seed for Batch Correction


#Combine SE adn ZE for later use in Batch correction application




#Tissue Specificity Section


#Gopher Enrichment


