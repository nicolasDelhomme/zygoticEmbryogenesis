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
#suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(UpSetR))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(SuperExactTest))

#' Function for selecting significant genes (FDR < 0.01, |log2FC| => 0.5)
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
}

#' * Graphics
# pal <- brewer.pal(12,"Paired")
# hpal <- colorRampPalette(c("blue","white","red"))(100)
# mar <- par("mar")

#' Import data  
#' * DEGs in ZE
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_UpAndDown_one-percent-FDR-05-log2fc-cutoff.rda"))
#' * DEGs in SE
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")
dds_SE <- res_sig_genes
#' * Expression of genes in SE
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda")  
# this dataset includes also TEs, but they are not many (there are 291 TEs and ~66000 genes)

#' # Plots   
#' 
#' ##' Number of expressed genes in SE 
#' Check number of expressed genes per time point also for the SE experiment- if there would be too much variation that would break the assumptions of DE analysis.  
#' 
# expressed genes per time point (check average expression of all the replicates)
expressedGenesSE_byStage <- sapply(split.data.frame(t(assay(dds)),colData(dds)$Stages), 
                                   function(x){colMeans(x)>0})
expressedGenesSE_byStage_nr <- apply(expressedGenesSE_byStage, 2, sum)

# genes with vst value >1 in any of the time points
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

expressedGenesSE_vst1_byStage <- sapply(split.data.frame(t(vst), colData(dds)$Stages), 
                                        function(x){colMeans(x)>1})
expressedGenesSE_vst1_byStage_nr <- apply(expressedGenesSE_vst1_byStage, 2, sum)

expressedGenesSE_vst5_byStage <- sapply(split.data.frame(t(vst), colData(dds)$Stages), 
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


#' ### Intersections of DEGs   

#' #### ZE vs FMG in each time point  
#' All DEGs genes  
# venn diagram is not very useful to compare that many datasets, therefore I will use upset plot
venn(lapply(genes_ZEFMG_all, function(x) rownames(x)), zcolor = "style", box = FALSE)
#' ZE vs FMG: all DEGs
upset(fromList(lapply(genes_ZEFMG_all, function(x) rownames(x))), sets = names(genes_ZEFMG_all), keep.order = TRUE, sets.bar.color = "grey")


#' ZE vs FMG: down-regulated genes
upset(fromList(lapply(genes_ZEFMG_down, function(x) rownames(x))), sets = names(genes_ZEFMG_down), keep.order = TRUE, sets.bar.color = "grey")
#' ZE vs FMG: up-regulated genes
upset(fromList(lapply(genes_ZEFMG_up, function(x) rownames(x))), sets = names(genes_ZEFMG_up), keep.order = TRUE, sets.bar.color = "grey")


########## How does union of the genes above (from every timepoint comparison) compares to Tissue_ZE_vs_FMG? CHECK

#' #### Venn diagram of DEGs at any point in SE, ZE,FMG and S (make union of all DEGs in a tissue regardless of time)  
# Pool together DEGs from all timepoint comparissons in each tissue
allDEGs_S <- unique(unlist(lapply(genes_S_down, rownames), use.names = FALSE), unlist(lapply(genes_S_up, rownames), use.names = FALSE))
allDEGs_FMG <- unique(unlist(lapply(genes_FMG_down, rownames), use.names = FALSE), unlist(lapply(genes_FMG_up, rownames), use.names = FALSE))
allDEGs_ZE <- unique(unlist(lapply(genes_ZE_down, rownames), use.names = FALSE), unlist(lapply(genes_ZE_up, rownames), use.names = FALSE))

# DEGs from SE
allDEGs_SE <- unique(unlist(lapply(dds_SE, rownames), use.names = FALSE))

#' Intersections of all DEGs regardless of time in different datasets/tissues
venn(list(seed = allDEGs_S, zygotic_embryo = allDEGs_ZE, FMG = allDEGs_FMG, somatic_embryo = allDEGs_SE), zcolor = "style", box = FALSE)
upset(fromList(list(seed = allDEGs_S, zygotic_embryo = allDEGs_ZE, FMG = allDEGs_FMG, somatic_embryo = allDEGs_SE)), keep.order = TRUE, sets.bar.color = "grey")

#' Can we get intersections of DEGs between all the stage comparisons of different kinds of the tissue?  
# combine lists of DEGs from all the tissues
DEGsListAll <- lapply(list(genes_S_all, genes_FMG_all, genes_ZE_all, dds_SE), function(x){
  lapply(x, rownames)})
DEGsListAll <- unlist(DEGsListAll, recursive = FALSE)
names(DEGsListAll) <- sub(pattern = "res_", replacement = "SE_", names(DEGsListAll))
save(DEGsListAll, file=here("analysis/DE/ZE-FMG-allStages_duplSsamples/DEGsListAll_ZEandSE.rda"))

# plot intersections arranged by number of DEGs in the intersection, decreasing
upset(fromList(DEGsListAll),
      sets.bar.color = "grey", 
      nsets = 21, 
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison", 
      mb.ratio = c(0.5, 0.5),
      nintersects = 80)

#' all datasets: arranged by size of the intersection, decreasing
upset(fromList(DEGsListAll),
      order.by = "freq",
      sets = names(DEGsListAll)[order(names(DEGsListAll))],
      keep.order = TRUE,
      sets.bar.color = "grey",
     # nsets = 21,
      mainbar.y.label = "Number of DEGs in the intersection",
      sets.x.label = "Number of DEGs in the comparison",
      mb.ratio = c(0.5, 0.5),
      nintersects = 80)



#' #### DEGS in ZE and FMG (separately for each tissue), through time

#' FMG: all DE genes
upset(fromList(lapply(genes_FMG_all, rownames)), keep.order = TRUE, sets.bar.color = "grey")
#' ZE: all DE genes
upset(fromList(lapply(genes_ZE_all, rownames)), keep.order = TRUE, sets.bar.color = "grey")

#' FMG: up-regulated genes
upset(fromList(lapply(genes_FMG_up, rownames)), keep.order = TRUE, sets.bar.color = "grey")
#' FMF: down-regulated genes
upset(fromList(lapply(genes_FMG_down, rownames)), keep.order = TRUE, sets.bar.color = "grey")

#' ZE: up-regulated genes
upset(fromList(lapply(genes_ZE_up, rownames)), keep.order = TRUE, sets.bar.color = "grey")
#' ZE: down-regulated genes
upset(fromList(lapply(genes_ZE_down, rownames)), keep.order = TRUE, sets.bar.color = "grey")


#' # Probability of multi-set intersections  
#' x = a list of character vectors, n = background population size
#' from the vignette: total = 18196; cis-eQTL gene sets were independently and randomly sampled from the population of 18,196 unique genes profiled in the eQTL study  
#' What is total? All genes in the spruce genome? Or only expressed? But in which experiment (SE/ZE)? What if you are comparing genes from two different genomes?
#' Should that be union of all unique expressed/DE gene names then?
total = 66000
resEST_allStageComp <- supertest(x = DEGsListAll, n = total)
save(resEST_allStageComp, file = here("analysis/DE/ZE-FMG-allStages_duplSsamples/resEST_allStageComp.rda"))
plot(resEST, sort.by = "size", degree = 2:15, keep.empty.intersections=FALSE, min.intersection.size=10)

#' # Expression of common genes in different tissues  
#' Gene DE between ZE and FMG at any time
# all_DEGs_ZEFMG <- unique(unlist(lapply(genes_ZEFMG_all, rownames), use.names = FALSE))
# 
# vsda <- varianceStabilizingTransformation(dds, blind = FALSE)
# vsta <- assay(vsda)
# vsta <- vsta - min(vsta)
# 
# vsta_ZEFMG <- vsta[all_DEGs_ZEFMG[all_DEGs_ZEFMG %in% rownames(vsta)], ]
# 
# heatmap.2(vsta_ZEFMG,labCol = all.colData$Stages, labRow = FALSE, trace = "none")  


#' #' 
#' #' ```{r empty,eval=FALSE,echo=FALSE}
#' #' ```
#' #' 
#' #' # Session Info
#' #' ```{r session info, echo=FALSE}
#' #' sessionInfo()
#' #' ```
#' 
