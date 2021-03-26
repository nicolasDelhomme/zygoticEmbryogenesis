#' ---
#' title: "Exploring DEGs in ZE, FMG and SE samples"
#' author: "Katja Stojkovic"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' 
#' # Setup  
#' 
#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(SuperExactTest))

#' Get helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))
source("~/Git/zygoticEmbryogenesis/Rtoolbox/src/plotEnrichedTreemap.R")

#' * Graphics
pal <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' Import data  
#' * DEGs in ZE
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_UpAndDown_one-percent-FDR-05-log2fc-cutoff.rda"))
#' * DEGs in SE
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")
ddsDEGs_SE <- res_sig_genes
rm(res_sig_genes)
#' * Expression of genes in SE
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda") 
dds_SE <- dds
rm(dds)

#' * Expression of genes in ZE
load(here("analysis/salmon/ZE-allStages_duplS_dds.rda"))
dds_ZE <- dds
rm(dds)

#' # Plots   
#' 
#' ##' Number of expressed genes in SE 
#' Check number of expressed genes per time point for the SE experiment- if there would be too much variation that would break the assumptions of DE analysis.  
#' 
# this dataset includes also TEs, but they are not many (there are 291 TEs and ~66000 genes)  
# exclude them:-----------------------------------????

# expressed genes per time point (check average expression of all the replicates)
expressedGenesSE_byStage <- sapply(split.data.frame(t(assay(dds_SE)),colData(dds_SE)$Stages), 
                                   function(x){colMeans(x)>0})
expressedGenesSE_byStage_nr <- apply(expressedGenesSE_byStage, 2, sum)

# genes with vst value >1 in any of the time points
vsd_SE <- varianceStabilizingTransformation(dds_SE, blind = TRUE)
vst_SE <- assay(vsd_SE)
vst_SE <- vst_SE - min(vst_SE)

expressedGenesSE_vst1_byStage <- sapply(split.data.frame(t(vst_SE), colData(dds_SE)$Stages), 
                                        function(x){colMeans(x)>1})
expressedGenesSE_vst1_byStage_nr <- apply(expressedGenesSE_vst1_byStage, 2, sum)

expressedGenesSE_vst5_byStage <- sapply(split.data.frame(t(vst_SE), colData(dds_SE)$Stages), 
                                        function(x){colMeans(x)>5})
expressedGenesSE_vst5_byStage_nr <- apply(expressedGenesSE_vst5_byStage, 2, sum)

# Plot number of genes expressed or having vst value >1 in different time points
barplot(rbind(expressedGenesSE_byStage_nr, expressedGenesSE_vst1_byStage_nr, expressedGenesSE_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        main = "Genes expressed at different values in SE", 
        sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"))
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")

#' ##' Number of expressed genes in ZE 
#' Check number of expressed genes per time point for the ZE experiment- if there would be too much variation that would break the assumptions of DE analysis.  
# 
# expressed genes per time point (check average expression of all the replicates)
expressedGenesZE_byStage <- sapply(split.data.frame(t(assay(dds_ZE)),colData(dds_ZE)$Time), 
                                 function(x){colMeans(x)>0})
expressedGenesZE_byStage_nr <- apply(expressedGenesZE_byStage, 2, sum)

# genes (average of replicates) with vst value >1 in any of the time points
vsd_ZE <- varianceStabilizingTransformation(dds_ZE, blind = TRUE)
vst_ZE <- assay(vsd_ZE)
vst_ZE <- vst_ZE - min(vst_ZE)

expressedGenesZE_vst1_byStage <- sapply(split.data.frame(t(vst_ZE), colData(dds_ZE)$Time), 
                                      function(x){colMeans(x)>1})
expressedGenesZE_vst1_byStage_nr <- apply(expressedGenesZE_vst1_byStage, 2, sum)

expressedGenesZE_vst5_byStage <- sapply(split.data.frame(t(vst_ZE), colData(dds_ZE)$Time), 
                                      function(x){colMeans(x)>5})
expressedGenesZE_vst5_byStage_nr <- apply(expressedGenesZE_vst5_byStage, 2, sum)

# Plot number of genes expressed or having vst value >1 in different time points
barplot(rbind(expressedGenesZE_byStage_nr, expressedGenesZE_vst1_byStage_nr, expressedGenesZE_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        main = "Genes expressed at different values in ZE", 
        sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"))
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")

#' ## Intersections of DEGs (upset plot)  
#' Set the commonly used parameters
upset_default <- function(x){upset(fromList(lapply(x, rownames)), sets = names(x), keep.order = TRUE, sets.bar.color = "grey", nsets = length(x), nintersects = NA)}
# nsets = "nr of sets", others only 5 will be shown,  
# nintersects = NA to display all of them, otherwise max 40 will be shown  
# for keep.order to work, sets need to be specified; an order of sets provided to sets is kept  

#' ### DEGS in ZE and FMG (separately for each tissue), through time  

#' FMG: all DE genes
upset_default(genes_FMG_all)
#' ZE: all DE genes
upset_default(genes_ZE_all)

#' FMG: up-regulated genes
upset_default(genes_FMG_up)
#' FMF: down-regulated genes
upset_default(genes_FMG_down)

#' ZE: up-regulated genes
upset_default(genes_ZE_up)
#' ZE: down-regulated genes
upset_default(genes_ZE_down)


#' ### ZE vs FMG in each time point  

#' ZE vs FMG: all DEGs
upset_default(genes_ZEFMG_all)
#' ZE vs FMG: down-regulated genes
upset_default(genes_ZEFMG_down)
#' ZE vs FMG: up-regulated genes
upset_default(genes_ZEFMG_up)


########## How does union of the genes above (from every timepoint comparison) compares to Tissue_ZE_vs_FMG? CHECK

#' ### Comparisson of DEGs at any point in SE, ZE,FMG and S (make union of all DEGs in a tissue regardless of time)  
# Pool together DEGs from all timepoint comparissons in each tissue
allDEGs_S <- unique(unlist(lapply(genes_S_all, rownames), use.names = FALSE))
allDEGs_FMG <- unique(unlist(lapply(genes_FMG_all, rownames), use.names = FALSE))
allDEGs_ZE <- unique(unlist(lapply(genes_ZE_all, rownames), use.names = FALSE))

# DEGs from SE
allDEGs_SE <- unique(unlist(lapply(ddsDEGs_SE, rownames), use.names = FALSE))

#' Intersections of all DEGs regardless of time in different datasets/tissues
DEGsListAllGeneral <- list(seed = allDEGs_S, zygotic_embryo = allDEGs_ZE, FMG = allDEGs_FMG, somatic_embryo = allDEGs_SE)
venn::venn(DEGsListAllGeneral, zcolor = "style", box = FALSE)
upset(fromList(DEGsListAllGeneral), sets = names(DEGsListAllGeneral), keep.order = TRUE, sets.bar.color = "grey", nsets = length(DEGsListAllGeneral), nintersects = NA)

#' Can we get intersections of DEGs between all the stage comparisons of different kinds of the tissue?  
# combine lists of DEGs from all the tissues
DEGsListAllComp <- lapply(list(genes_S_all, genes_FMG_all, genes_ZE_all, ddsDEGs_SE), function(x){
  lapply(x, rownames)})
DEGsListAllComp <- unlist(DEGsListAllComp, recursive = FALSE)
names(DEGsListAllComp) <- sub(pattern = "res_", replacement = "SE_", names(DEGsListAllComp))
save(DEGsListAllComp, file=here("analysis/DE/ZE-FMG-allStages_duplSsamples/DEGsListAllComp_ZEandSE.rda"))

# plot intersections arranged by number of DEGs in the intersection, decreasing
upset(fromList(DEGsListAllComp),
      sets = names(DEGsListAllComp),
      keep.order = TRUE,
      order.by = "freq",
      sets.bar.color = "grey", 
      nsets = length(DEGsListAllComp), 
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison", 
      mb.ratio = c(0.5, 0.5),
      nintersects = 80,
      number.angles = 45)


#' ## Intersections of DEGs with probability (SuperExactTest)  
#' How likely is that intersection of this size would appear by chance (because of the size of the group)? 
#' Is it expected that this many genes would appear the same, comparing two groups of genes with certain background population size?  
#' Calculated p value tells us one-tail probability of observing equal to or larger than the number of intersect items.
#' 
#' How to choose a background population size (n)?  
#' From the vignette: in their case n = total = 18196; cis-eQTL gene sets were independently and randomly sampled from the population of 18,196 unique genes profiled in the eQTL study  
#' 
#' We should use a number of expressed genes in the dataset (experiment). When comparing two experiments (SE & ZE) use a union of expressed genes.  
#' Comment (Nico): It is easier to assume that the population consists of any gene expressed in that dataset (as we do for GO). It might vary slightly between tissues and time point, but not drastically as shown above. 
#' So it won’t affect the stats much.  
  
#' Calculate number of expressed genes in SE and ZE:  
# i) SE  
# all the genes with expression > 0 in at least one time point
expressedGenesSE_anyStage <- apply(expressedGenesSE_byStage, 1, any)
# how many are they?
SEbackground_0 <- sum(expressedGenesSE_anyStage)

# all the genes with vst > 01 in at least one time point
expressedGenesSE_vst1_anyStage <- apply(expressedGenesSE_vst1_byStage, 1, any)
# how many are they?
SEbackground_1 <- sum(expressedGenesSE_vst1_anyStage)

# ii) ZE  
# all the genes with expression > 0 in at least one time point
expressedGenesZE_anyStage <- apply(expressedGenesZE_byStage, 1, any)
# how many are they?
ZEbackground_0 <- sum(expressedGenesZE_anyStage)

# all the genes with vst > 1 in at least one time point
expressedGenesZE_vst1_anyStage <- apply(expressedGenesZE_vst1_byStage, 1, any)
# how many are they?
ZEbackground_1 <- sum(expressedGenesZE_vst1_anyStage)

#' Number of expressed genes in SE and ZE experiments (raw expression > 0 or vst >1; average of replicates in at least one time)
pander(data.frame("background_0" = c(ZEbackground_0, SEbackground_0), "background_1" = c(ZEbackground_1, SEbackground_1), row.names = c("ZE", "SE")))

#' VST expression is relative to the dataset/experiment, as we subtract the minimium to get vst values and it is therefore not directly comparable between experiments. 
#' Could we use a different measure of expression?  
#' * Estimate noise and filter expressed genes (check the script in UPSCb-comon)  
#' * Count all the genes taken for doing DE analysis; that would be all the genes not discarded (padj is not NA) by independent filtering in DESeq2.  
#'   
#' Import not filtered DE results for SE and ZE in order to get number of expressed genes from DESeq objects  
#' ZE
load(here("data/analysis/DE/ZE-FMG-allStages_duplSsamples/notfilteredDEGs.rda"))
# Number of genes (whose padj is not NA) reported in results of each of the comparison is not the same
listZEstages <- lapply(list(resIntercept, B2.S_B1.S, B3.S_B2.S, 
                     B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG,
                     B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE),
                     function(x) {x[!is.na(x$padj), ]})
names(listZEstages) <- c("Intercept", "B2.S_B1.S", "B3.S_B2.S", 
                         "B4.FMG_B3.S", "B5.FMG_B4.FMG", "B6.FMG_B5.FMG", "B7.FMG_B6.FMG", "B8.FMG_B7.FMG", "B9.FMG_B8.FMG",
                         "B4.ZE_B3.S", "B5.ZE_B4.ZE", "B6.ZE_B5.ZE", "B7.ZE_B6.ZE", "B8.ZE_B7.ZE", "B9.ZE_B8.ZE")

NrExpGeneZEstage <- sapply(listZEstages, nrow)

pander(NrExpGeneZEstage)
pander(table(NrExpGeneZEstage))
upset_default(listZEstages)

# Names of all the expressed genes in ZE (union of expressed genes in all the comparisons):
NameExpGeneZE <- unique(unlist(lapply(listZEstages, rownames), use.names = FALSE))
# Number of expressed genes in zygotic embryogenesis (after independent filtering) is therefore:
NrExpGeneZE <- length(NameExpGeneZE)
NrExpGeneZE
# This represents % of all genes (note: number of genes before filtering is the same in any comparison: intercept, B2.S_B3.S, ...)
NrExpGeneZE/nrow(resIntercept)*100

# Could we use intercept as the approximation of all the expressed genes in the experiment? Check also in other experiments!
nrow(listZEstages$Intercept)
# check the threshold for independent filtering: % of removed genes and the cutoff value for expression (baseMean)
metadata(resIntercept)$filterThreshold

# remove unnecessary objects
rm(resIntercept, 
   B2.S_B1.S, B3.S_B2.S, 
   B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG, 
   B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE,
   B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG)


#' SE  
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered.rda")
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered_onlyIntercept.rda")

listSEstages <- res_list
listSEstages[["Intercept"]] <- resInterceptSE

# exclude TEs
listSEstages <- lapply(listSEstages, function(x){
  x[grepl("^MA_", rownames(x)), ]
})

# filter for expressed genes
listSEstages <- lapply(listSEstages, function(x){
  x[!is.na(x$padj), ]
})

NrExpGeneSEstage <- sapply(listSEstages, nrow)

pander(NrExpGeneSEstage)
pander(table(NrExpGeneSEstage))
upset_default(listSEstages)

# Names of all the expressed genes in SE (union of expressed genes in all the comparisons):
NameExpGeneSE <- unique(unlist(lapply(listSEstages, rownames), use.names = FALSE))
# Number of expressed genes in zygotic embryogenesis (after independent filtering) is therefore:
NrExpGeneSE <- length(NameExpGeneSE)
NrExpGeneSE
# This represents % of all genes (note: number of genes before filtering is the same in any comparison: intercept, S2 vs S1, ... = 66069)
NrExpGeneSE/66069*100

# Could we use intercept as the approximation of all the expressed genes in the experiment? Check also in other experiments!
nrow(listSEstages$Intercept)
# check the threshold for independent filtering: % of removed genes and the cutoff value for expression (baseMean)
metadata(listSEstages$Intercept)$filterThreshold

# remove unnecessary objects
rm(resInterceptSE, res_list)


# SuperExactTest plots----------------------------------------------------  


#' Compare setup & SET plots:  
#' * "intersection" of one set differs: in setup plot only specific genes for that group are represented, but in SET plot number of all the genes in the group is represented  
#' * there are some minor differences in number of genes in the intersections - how is this possible?  

#' ### DEGS in ZE and FMG (separately for each tissue), through time  
#' Define size of background population
n = NrExpGeneZE

#' FMG: all DE genes
setFMG_all <- supertest(lapply(genes_FMG_all, rownames), n = n)
plot(setFMG_all, 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape", 
     bar.split= c(800, 3800), color.on = "black", color.off = "white")
dev.off()
# plot P values
# plot(setFMG_all$P.value)
# abline(h=0.05, col = "red")
# summary(setFMG_all)$Table[summary(setFMG_all)$P.value < 0.05, ]

#' ZE: all DE genes
setZE_all <- supertest(lapply(genes_ZE_all, rownames), n = n)
dev.new()
plot(setZE_all, 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape", 
     bar.split= c(3000, 6000), color.on = "black", color.off = "white")
dev.off()

# plot P values
# plot(setZE_all$P.value)
# abline(h=0.05, col = "red")
# summary(setZE_all)$Table[summary(setZE_all)$P.value < 0.05, ]

#' FMG: up-regulated genes
dev.new()
plot(supertest(lapply(genes_FMG_up, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape", 
     bar.split= c(600, 1400), color.on = "black", color.off = "white")
dev.off()
#' FMG: down-regulated genes
dev.new()
plot(supertest(lapply(genes_FMG_down, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     bar.split= c(500, 2300), color.on = "black", color.off = "white")
dev.off()

#' ZE: up-regulated genes
dev.new()
plot(supertest(lapply(genes_ZE_up, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     bar.split= c(1000, 3000),color.on = "black", color.off = "white")
dev.off

#' ZE: down-regulated genes
dev.new()
plot(supertest(lapply(genes_ZE_down, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     bar.split= c(1500, 3000), color.on = "black", color.off = "white")
dev.off()
     
     
#' ### ZE vs FMG in each time point  

#' ZE vs FMG: all DEGs
plot(supertest(lapply(genes_ZEFMG_all, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     color.on = "black", color.off = "white")
#' ZE vs FMG: down-regulated genes
plot(supertest(lapply(genes_ZEFMG_down, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     color.on = "black", color.off = "white")
#' ZE vs FMG: up-regulated genes
plot(supertest(lapply(genes_ZEFMG_up, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     color.on = "black", color.off = "white")

#' ### Comparisson of DEGs at any point in SE, ZE,FMG and S (make union of all DEGs in a tissue regardless of time)  
#' Check differences in the background population (all the expressed genes)  
#' #### raw and vst counts  
# exp > 0
expressedGenesZE_anyStage_names <- names(expressedGenesZE_anyStage)[expressedGenesZE_anyStage]
expressedGenesSE_anyStage_names <- names(expressedGenesSE_anyStage)[expressedGenesSE_anyStage]
# vst > 1
expressedGenesZE_vst1_anyStage_names <- names(expressedGenesZE_vst1_anyStage)[expressedGenesZE_vst1_anyStage]
expressedGenesSE_vst1_anyStage_names <- names(expressedGenesSE_vst1_anyStage)[expressedGenesSE_vst1_anyStage]
# compare genes expressed at different level
venn::venn(list(SE_0 = expressedGenesSE_anyStage_names, ZE_0 = expressedGenesZE_anyStage_names, 
                SE_1 = expressedGenesSE_vst1_anyStage_names, ZE_1 = expressedGenesZE_vst1_anyStage_names), 
           zcolor = "style", box = FALSE)

#' Use union of genes expressed in SE and ZE projects when comparing groups of DEGs between them.  
# exp > 0
expressedGenesZESE_anyStage_names <- unique(c(expressedGenesZE_anyStage_names, expressedGenesSE_anyStage_names))
# vst > 1
expressedGenesZESE_vst1_anyStage_names <- unique(c(expressedGenesZE_vst1_anyStage_names, expressedGenesSE_vst1_anyStage_names))
# compare
venn::venn(list(ZESE_0 = expressedGenesZESE_anyStage_names, ZESE_1 = expressedGenesZESE_vst1_anyStage_names),
           zcolor = "style", box = FALSE)

#' #### Nr of genes after independent filtering in DESeq2  
venn::venn(list(SE = NameExpGeneSE, ZE = NameExpGeneZE),
           zcolor = "style", box = FALSE)

#' Union of genes expressed in SE and ZE
NameExpGeneSEZE <- union(NameExpGeneSE, NameExpGeneZE)
#' Number of genes in the union
length(NameExpGeneSEZE)

#' Intersections of all DEGs regardless of time in different datasets/tissues  
#' Define size of background population (union of unique expressed genes in SE and ZE)
n = length(NameExpGeneSEZE)

plot(supertest(DEGsListAllGeneral, n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
     color.on = "black", color.off = "white")

#' Can we get intersections of DEGs between all the stage comparisons of different kinds of the tissue? - TAKES LONG TIME
# SET_DEGsListAllComp <- supertest(DEGsListAllComp, n = n)
# save(SET_DEGsListAllComp, file = here("analysis/DE/ZE-FMG-allStages_duplSsamples/SET_allStageComp.rda"))
# 
# # plot by degree, degree = 1-3
# plot(SET_DEGsListAllComp,
#      sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
#      color.on = "black", color.off = "white", degree = c(1:3), min.intersection.size = 10, sort.by = c("degree", "size"))

# # plot by degree, degree = 4-10
# plot(SET_DEGsListAllComp,
#      sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
#      color.on = "black", color.off = "white", degree = c(4:10), min.intersection.size = 10, sort.by = c("degree", "size"))
# 
# # plot by degree, degree = 11-21
# plot(SET_DEGsListAllComp,
#      sort.by = "size", keep.empty.intersections = FALSE, Layout = "landscape",
#      color.on = "black", color.off = "white", degree = c(11:21), min.intersection.size = 10, sort.by = c("degree", "size"))

# Heatmaps and treemaps--------------------------------------------

#' # Expression of common genes in different tissues  
#' ## Design aware vst  
#' SE  
# vsda_SE <- varianceStabilizingTransformation(dds_SE, blind = FALSE)
# vsta_SE <- assay(vsda_SE)
# vsta_SE <- vsta_SE - min(vsta_SE)
# write.csv(vsta_SE, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", quote = FALSE)
vsta_SE <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", header = TRUE, row.names = 1)

#' ZE  
# vsda_ZE <- varianceStabilizingTransformation(dds_ZE, blind = FALSE)
# vsta_ZE <- assay(vsda_ZE)
# vsta_ZE <- vsta_ZE - min(vsta_ZE)
# write.csv(vsta_ZE, file=here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), quote = FALSE)
vsta_ZE <- read.csv(here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), header = TRUE, row.names = 1)

#' Scale the data
vstaSE_scaled <- t(scale(t(vsta_SE)))
vstaZE_scaled <- t(scale(t(vsta_ZE)))

#' ## Expression of DEGs  
#' 
#' ### DEGs in ZE, but not in SE  
DEGs_ZE_notSE <- setdiff(union(allDEGs_S, allDEGs_ZE, allDEGs_FMG), allDEGs_SE)

# without scaling
# heatmap.2(vsta_ZE[DEGs_ZE_notSE, ], labCol = colData(dds_ZE)$Replicate, labRow = FALSE, trace = "none", main = "ZE vs SE")
# scaled data
# heatmap.2(vsta_ZE[DEGs_ZE_notSE, ], scale = "row", labCol = colData(dds_ZE)$Replicate, labRow = FALSE, trace = "none", main = "ZE vs SE, scaled")

pdf(file = here("analysis/figures/DEGs_ZE_notSE_Replicate.pdf"), width = 15, height = 10)
#png(file = here("analysis/figures/DEGs_ZE_notSE.png"), width = 1000, height = 800, pointsize = 23)
#par(mar = c(20, 10, 10, 5)) # why it does not change the margins??
heatmap.2(vstaZE_scaled[DEGs_ZE_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs unique to zygotic embryogenesis", 
          xlab = "time points", 
          ylab = paste(length(DEGs_ZE_notSE), "genes"))
dev.off()

#' How are these genes expressed in SE?  
#' * How to deal with genes not expressed in SE and therefore NaN?  
#' * Keep the order of the genes the same when comparing their expression in two datasets/tissues  
#' 

#' Functional enrichment  
enrDEGs_ZE_notSE <- gopher(sub(".1$", "", DEGs_ZE_notSE),
                           task = list("go","mapman"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_ZE_notSE, file = here("analysis/gopher/enrDEGs_ZE_notSE.rda"))

#' ### DEGs in ZEm (B4-B9), but not in SE
DEGs_ZEm_notSE <- setdiff(allDEGs_ZE, allDEGs_SE)

pdf(file = here("analysis/figures/DEGs_ZEm_notSE_Replicate.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_ZEm_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(length(DEGs_ZEm_notSE), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_ZEm_notSE <- gopher(sub(".1$", "", DEGs_ZEm_notSE),
                           task = list("go","mapman"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_ZEm_notSE, file = here("analysis/gopher/enrDEGs_ZEm_notSE.rda"))

#' ### DEGs in FMG (B4-B9), but not in SE
DEGs_FMG_notSE <- setdiff(allDEGs_FMG, allDEGs_SE)

pdf(file = here("analysis/figures/DEGs_FMG_notSE_Replicate.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_FMG_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(length(DEGs_FMG_notSE), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_FMG_notSE <- gopher(sub(".1$", "", DEGs_FMG_notSE),
                            task = list("go","mapman"), 
                            background = sub(".1$", "", NameExpGeneSEZE), 
                            url="pabies")
save(enrDEGs_FMG_notSE, file = here("analysis/gopher/enrDEGs_FMG_notSE.rda"))

#' ### DEGs in SE, but not in ZE
DEGs_SE_notZE <- setdiff(allDEGs_SE, union(allDEGs_S, allDEGs_ZE, allDEGs_FMG))

pdf(file = here("analysis/figures/DEGs_SE_notZE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[DEGs_SE_notZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to somatic embryogenesis", 
          xlab = "time points",
          ylab = paste(length(DEGs_SE_notZE), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_SE_notZE <- gopher(sub(".1$", "", DEGs_SE_notZE),
                          task = list("go","mapman"), 
                          background = sub(".1$", "", NameExpGeneSEZE), 
                          url="pabies")
save(enrDEGs_SE_notZE, file = here("analysis/gopher/enrDEGs_SE_notZE.rda"))

#' ### DEGs in SE and ZE
DEGs_SEandZE <- intersect(allDEGs_SE, union(allDEGs_S, allDEGs_ZE, allDEGs_FMG))

#' Plot expression of the same genes in SE and ZE dataset.  
# check if genes are in the same order in both datasets
all(rownames(vstaZE_scaled[DEGs_SEandZE, ]) == rownames(vstaSE_scaled[DEGs_SEandZE, ]))
# join them by cbind ### how do you then set other parameters for heatmap? e.g. clustering of the samples inside its own dataset?  
# pdf(file = here("analysis/figures/DEGs_SEandZE.pdf"), width = 15, height = 10)
# heatmap.2(cbind(vstaSE_scaled[DEGs_SEandZE, ], vstaZE_scaled[DEGs_SEandZE, ]), 
#           trace = "none",
#           distfun = function(x) as.dist(1-cor(t(x))),
#           hclustfun = function(X){hclust(X,method="ward.D")},
#           #ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
#           col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
#           labCol = c(colData(dds_SE)$Stages, colData(dds_ZE)$Replicate),
#           labRow = NA, 
#           main = "DEGs common to SE and ZE", 
#           xlab = "time points",
#           ylab = paste(length(DEGs_SEandZE), "genes"))
# dev.off()

# make two plots and compare them side by side  
pdf(file = here("analysis/figures/DEGs_SEandZE_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[DEGs_SEandZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages,
          labRow = NA, 
          main = "DEGs common to SE and ZE (expression in SE)", 
          xlab = "time points",
          ylab = paste(length(DEGs_SEandZE), "genes"))
dev.off()

pdf(file = here("analysis/figures/DEGs_SEandZE_inZE.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_SEandZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs common to SE and ZE (expresion in ZE)", 
          xlab = "time points", 
          ylab = paste(length(DEGs_SEandZE), "genes"))
dev.off()

# do we want to fix order of the samples in this case to compare expression profile of the genes by the timeline in two datasets?  

#' Functional enrichment  
enrDEGs_SEandZE <- gopher(sub(".1$", "", DEGs_SEandZE),
                           task = list("go","mapman"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_SEandZE, file = here("analysis/gopher/enrDEGs_SEandZE.rda"))


#' ### DEGs in ZE compared to FMG at each time point  
#' make union of DEGs expressed at each time point
all_DEGs_ZEFMG <- unique(unlist(lapply(genes_ZEFMG_all, rownames), use.names = FALSE))

pdf(file = here("analysis/figures/DEGs_ZEvsFMG.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[all_DEGs_ZEFMG, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "Genes DE in ZE compared to FMG at any time point", 
          xlab = "time points", 
          ylab = paste(length(all_DEGs_ZEFMG), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_allZEFMG <- gopher(sub(".1$", "", all_DEGs_ZEFMG),
                        task = list("go","mapman"), 
                        background = sub(".1$", "", NameExpGeneZE), 
                        url="pabies")
save(enrDEGs_allZEFMG, file = here("analysis/gopher/enrDEGs_allZEFMG.rda"))

#' ### Genes DE in B7 compared to B6  
# make a union of genes DEGs between B7 and B6 in ZE or FMG
DEGs_B7vsB6 <- unique(c(rownames(genes_FMG_all$B7.FMG_B6.FMG), rownames(genes_ZE_all$B7.ZE_B6.ZE)))
#' How many DEGs do ZE and FMG have in common in this time point comparison?
dev.new()
venn::venn(list(rownames(genes_FMG_all$B7.FMG_B6.FMG), rownames(genes_ZE_all$B7.ZE_B6.ZE)), zcolor = "style", box = FALSE)

# divide genes in three groups  
B7vB6_FMGuniq <- setdiff(rownames(genes_FMG_all$B7.FMG_B6.FMG), rownames(genes_ZE_all$B7.ZE_B6.ZE))
B7vB6_ZEuniq <- setdiff(rownames(genes_ZE_all$B7.ZE_B6.ZE), rownames(genes_FMG_all$B7.FMG_B6.FMG))
B7vB6_FMGZEint <- intersect(rownames(genes_ZE_all$B7.ZE_B6.ZE), rownames(genes_FMG_all$B7.FMG_B6.FMG))

# prepare vector with names of the groups, reorder it
B7vB6_groupByTissue <- rep(c("FMG", "ZE", "FMG&ZE"), c(length(B7vB6_FMGuniq), length(B7vB6_ZEuniq), length(B7vB6_FMGZEint)))
names(B7vB6_groupByTissue) <- c(B7vB6_FMGuniq, B7vB6_ZEuniq, B7vB6_FMGZEint)
B7vB6_groupByTissue <- B7vB6_groupByTissue[DEGs_B7vsB6]

pdf(file = here("analysis/figures/DEGs_B7vsB6.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_B7vsB6, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          RowSideColors = c("blue", "yellow", "green")[as.integer(factor(B7vB6_groupByTissue))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "Genes DE in B7 compared to B6 in ZE or FMG", 
          xlab = "time points", 
          ylab = paste(length(DEGs_B7vsB6), "genes"))
# legend("bottomleft",      
#        legend = levels(factor(B7vB6_groupByTissue)),
#        col = c("blue", "green", "yellow")[1:length(levels(factor(B7vB6_groupByTissue)))], 
#        fill = c("blue", "green", "yellow")[1:length(levels(factor(B7vB6_groupByTissue)))])
dev.off()

#' How are these genes expressed in SE?  
#' How many of these genes is expressed in SE?  
sum(DEGs_B7vsB6 %in% NameExpGeneSE)
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% NameExpGeneSE)/length(DEGs_B7vsB6)*100
#' How many of these genes is DE in SE? 
sum(DEGs_B7vsB6 %in% allDEGs_SE) 
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% allDEGs_SE)/length(DEGs_B7vsB6)*100

# What about vst?
sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled))
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled))/length(DEGs_B7vsB6)*100

#' Plot expression of genes DE between B7 and B6 in SE - not the same order of genes!
# B7vsB6_inSE <- vstaSE_scaled[DEGs_B7vsB6, ]
# B7vsB6_inSE <- B7vsB6_inSE[!is.na(B7vsB6_inSE)]
# pdf(file = here("analysis/figures/DEGs_B7vsB6_inSE.pdf"), width = 15, height = 10)
# heatmap.2(vstaSE_scaled[DEGs_B7vsB6, ], ###################### where do NaN values come from?
#           trace = "none",
#           distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
#           hclustfun = function(X){hclust(X,method="ward.D")},
#           ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
#           RowSideColors = c("blue", "yellow", "green")[as.integer(factor(B7vB6_groupByTissue))],
#           col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
#           labCol = colData(dds_SE)$Stages, 
#           labRow = NA, 
#           main = "Genes DE in B7 compared to B6 in ZE or FMG (expression in SE)", 
#           xlab = "time points", 
#           ylab = paste(sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled)), "genes"))
# dev.off()

#' Functional enrichement  
enrDEGs_B7vB6 <- gopher(sub(".1$", "", DEGs_B7vsB6),
                        task = list("go","mapman"), 
                        background = sub(".1$", "", NameExpGeneZE), 
                        url="pabies")
save(enrDEGs_B7vB6, file = here("analysis/gopher/enrDEGs_B7vB6.rda"))

#' ### Plot treemaps  
treemapList <- list(enrDEGs_allZEFMG, enrDEGs_B7vB6, enrDEGs_SEandZE, enrDEGs_ZE_notSE, enrDEGs_ZEm_notSE, enrDEGs_FMG_notSE, enrDEGs_SE_notZE)
names(treemapList) <- c("enrDEGs_allZEFMG", "enrDEGs_B7vB6", "enrDEGs_SEandZE", "enrDEGs_ZE_notSE", "enrDEGs_ZEm_notSE", "enrDEGs_FMG_notSE", "enrDEGs_SE_notZE")

lapply(names(treemapList), function(x){
  pdf(file = paste0(here("analysis/gopher/"), x, "_treemap.pdf"), width = 10, height = 6)
  plotEnrichedTreemap(treemapList[[x]], enrichment = "go", namespace = "none")
  plotEnrichedTreemap(treemapList[[x]], enrichment = "mapman")
  dev.off() 
})


#' #' 
#' #' ```{r empty,eval=FALSE,echo=FALSE}
#' #' ```
#' #' 
#' #' # Session Info
#' #' ```{r session info, echo=FALSE}
#' #' sessionInfo()
#' #' ```
#' 
